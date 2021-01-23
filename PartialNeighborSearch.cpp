#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <chrono>
#include <fstream>

using namespace std;
using namespace chrono;

#define N 8
#define M (N*N + 1) //The number of different magnet levels for given N
#define seed 19950724
#define Inf numeric_limits<double>::infinity()
default_random_engine generator(seed);
uniform_real_distribution<double> uniformRealDistribution(0, 1);
uniform_int_distribution<int> uniformIntDistribution(0, N - 1);

class Ising {
public:
    int matrix[N][N];

    void Random() {
        for (auto & i : matrix) {
            for (int & j : i) {
                j = (uniformRealDistribution(generator) < 0.5) ? -1 : 1;
            }
        }
    }

    int Magnet() {
        int ans = 0;
        for (auto & i : matrix) {
            for (int j : i) {
                ans += j;
            }
        }
        return ans;
    }

    int Energy() {
        int ans = 0;
        for (int i = 0; i < N - 1; i++) {
            for (int j = 0; j < N - 1; j++) {

                ans = ans - matrix[i][j] * matrix[i + 1][j];
                ans = ans - matrix[i][j] * matrix[i][j + 1];
            }
        }
        for (int i = 0; i < N - 1; i++) {
            ans = ans - matrix[i][N - 1] * matrix[i + 1][N - 1];
        }
        for (int j = 0; j < N - 1; j++) {
            ans = ans - matrix[N - 1][j] * matrix[N - 1][j + 1];
        }
        return ans;
    }

    double Boltzmann(double temp) {
        double ans = exp(- Energy() / temp);
        return ans;
    }

    void Print(double temp = 1) {
        for (auto & i : matrix) {
            for (int j : i) {
                cout << j << " ";
            }
            cout << endl;
        }
        cout << "Magnet = " << Magnet();
        cout << endl;
        cout << "Energy = " << Energy();
        cout << endl;
        cout << "For temperature = " << temp << ", the Boltzmann = " << Boltzmann(temp);
        cout << endl;
    }

    void Flip(int x, int y) {
        matrix[x][y] = - matrix[x][y];
    }

    void UpdateConfiguration(double temp, int x = -1, int y = -1) {
        if (x == -1) {
            x = uniformIntDistribution(generator);
        }
        if(y == -1){
            y = uniformIntDistribution(generator);
        }
        double boltzmann_now = Boltzmann(temp);
        Flip(x, y);
        double boltzmann_new = Boltzmann(temp);
        Flip(x, y);
        double ratio = boltzmann_new / boltzmann_now;
        double Prob = (1 < ratio) ? 1 : ratio;
        double prob = uniformRealDistribution(generator);
        if (Prob > prob) {
            Flip(x, y);
        }
    }

    double UpdateProbability(double temp, int x, int y) {
        if (x == -1) {
            x = uniformIntDistribution(generator);
        }
        if(y == -1){
            y = uniformIntDistribution(generator);
        }
        double boltzmann_now = Boltzmann(temp);
        Flip(x, y);
        double boltzmann_new = Boltzmann(temp);
        Flip(x, y);
        double ratio = boltzmann_new / boltzmann_now;
        double Prob = (1 < ratio) ? 1 : ratio;
        return Prob;
    }
};

double *Magnet_probability(double temp) {
    auto *p_actual = new double[M];
    for (int i = 0; i < M; i++) {
        p_actual[i] = 0;
    }
    for (int i = 0; i < pow(2, (N * N)); i++) {
        int ij = i;
        Ising S{};
        for (int j = 0; j < (N * N); j++) {
            int x = j / N;
            int y = j % N;
            S.matrix[x][y] = (ij % 2) * 2 - 1;
            ij /= 2;
        }
        p_actual[(S.Magnet() + N * N) / 2] += S.Boltzmann(temp);
    }
    double p_sum = 0;
    for (int i = 0; i < M; i++) {
        p_sum += p_actual[i];
    }
    for (int i = 0; i < M; i++) {
        p_actual[i] /= p_sum;
    }
    return p_actual;
}

double *Metropolis(double temp, int Iter) {
    Ising Si{};
    Si.Random();
    double time_sum = 0;
    auto *p_sampled = new double[M];
    for (int i = 0; i < Iter; i++) {
        auto start1 = steady_clock::now();

        Si.UpdateConfiguration(temp);

        p_sampled[(Si.Magnet() + N * N) / 2] += 1;

        auto end1 = steady_clock::now();
        double time1 = duration_cast<nanoseconds>(end1 - start1).count();
        time_sum += time1;
    }
    double p_sum = 0;
    for (int i = 0; i < M; i++) {
        p_sum += p_sampled[i];
    }
    for (int i = 0; i < M; i++) {
        p_sampled[i] /= p_sum;
    }

    double *p_actual = Magnet_probability(temp);
    double tvd = 0;
    for (int i = 0; i < M; i++) {
        tvd += abs(p_sampled[i] - p_actual[i]) / 2;
    }
    auto *ans = new double[2];
    ans[0] = tvd;
    ans[1] = time_sum;
    return ans;
}

double *Metropolis_AlternatingChains(double temp, int Iter, int K) {
    Ising Si{};
    Si.Random();
    double time_sum = 0;
    auto *p_sampled = new double[M];
    for (int i = 0; i < (Iter / N / K); i++) {
        auto start1 = steady_clock::now();
        for(int j = 0 ; j < N; j ++) {
            for(int k = 0; k < K; k++){
                Si.UpdateConfiguration(temp, j);
                p_sampled[(Si.Magnet() + N * N) / 2] += 1;
            }
        }
        auto end1 = steady_clock::now();
        double time1 = duration_cast<nanoseconds>(end1 - start1).count();
        time_sum += time1;
    }
    double p_sum = 0;
    for (int i = 0; i < M; i++) {
        p_sum += p_sampled[i];
    }
    for (int i = 0; i < M; i++) {
        p_sampled[i] /= p_sum;
    }
    double *p_actual = Magnet_probability(temp);
    double tvd = 0;
    for (int i = 0; i < M; i++) {
        tvd += abs(p_sampled[i] - p_actual[i]) / 2;
    }
    auto *ans = new double[2];
    ans[0] = tvd;
    ans[1] = time_sum;
    return ans;
}

double *RejectionFree(double temp, int Iter) {
    Ising Si{};
    Si.Random();
    double time_sum = 0;
    auto *p_sampled = new double[M];
    for (int i = 0; i < Iter; i++) {
        double mean_prob = 0;
        int max_time = 0;
        double d_min = Inf;
        int xx, yy;
        for (int x = 0; x < N; x++) {
            for (int y = 0; y < N; y++) {
                auto start1 = steady_clock::now();
                double prob = Si.UpdateProbability(temp, x, y);
                mean_prob += prob / (N * N);
                double r = uniformRealDistribution(generator);
                double d = -log(r) / prob;
                if (d < d_min) {
                    d_min = d;
                    xx = x;
                    yy = y;
                }
                auto end1 = steady_clock::now();
                int time1 = duration_cast<nanoseconds>(end1 - start1).count();
                max_time = (max_time > time1) ? max_time : time1;
            }
        }
        time_sum += max_time;
        int mm = 1;
        double sum = 0;
        double rr = uniformRealDistribution(generator);
        for(; mm < 1e6; mm++)
        {
            sum += pow(1 - mean_prob, mm-1) * mean_prob;
            if(sum > rr)
            {
                break;
            }
        }
        auto start2 = steady_clock::now();
        p_sampled[(Si.Magnet() + N * N) / 2] += mm;
        Si.Flip(xx, yy);
        auto end2 = steady_clock::now();
        int time2 = duration_cast<nanoseconds>(end2 - start2).count();
        time_sum += time2;
    }
    double p_sum = 0;
    for (int i = 0; i < M; i++) {
        p_sum += p_sampled[i];
    }
    for (int i = 0; i < M; i++) {
        p_sampled[i] /= p_sum;
    }
    double tvd = 0;
    double *p_actual = Magnet_probability(temp);
    for (int i = 0; i < M; i++) {
        tvd += abs(p_sampled[i] - p_actual[i]) / 2;
    }
    auto *ans = new double[2];
    ans[0] = tvd;
    ans[1] = time_sum;
    return ans;
}

double *RejectionFree_AlternatingChains(double temp, int Iter, int K) {
    Ising Si{};
    Si.Random();
    double time_sum = 0;
    auto *p_sampled = new double[M];
    for (int i = 0; i < Iter / N; i++) {
        for(int j = 0; j < N ; j ++) {
            int k = 0;
            while(k < K) {
                double mean_prob = 0;
                int max_time = 0;
                double d_min = Inf;
                int yy;
                for (int y = 0; y < N; y++) {
                    auto start1 = steady_clock::now();
                    double prob = Si.UpdateProbability(temp, j, y);
                    mean_prob += prob / N;
                    double r = uniformRealDistribution(generator);
                    double d = -log(r) / prob;
                    if (d < d_min) {
                        d_min = d;
                        yy = y;
                    }
                    auto end1 = steady_clock::now();
                    int time1 = duration_cast<nanoseconds>(end1 - start1).count();
                    max_time = (max_time > time1) ? max_time : time1;
                }
                time_sum += max_time;
                int mm = 1;
                double sum = 0;
                double rr = uniformRealDistribution(generator);
                for(; mm < 1e6; mm++)
                {
                    sum += pow(1 - mean_prob, mm-1) * mean_prob;
                    if(sum > rr)
                    {
                        break;
                    }
                }

                auto start2 = steady_clock::now();
                if (k + mm  > K) {
                    p_sampled[(Si.Magnet() + N * N) / 2] += K - k;
                    k = K;
                }
                else {
                    p_sampled[(Si.Magnet() + N * N) / 2] += mm;
                    k += mm;
                    Si.Flip(j, yy);
                }
                auto end2 = steady_clock::now();
                int time2 = duration_cast<nanoseconds>(end2 - start2).count();
                time_sum += time2;
            }
        }
    }
    double p_sum = 0;
    for (int i = 0; i < M; i++) {
        p_sum += p_sampled[i];
    }
    for (int i = 0; i < M; i++) {
        p_sampled[i] /= p_sum;
    }
    double tvd = 0;
    double *p_actual = Magnet_probability(temp);
    for (int i = 0; i < M; i++) {
        tvd += abs(p_sampled[i] - p_actual[i]) / 2;
    }
    auto *ans = new double[2];
    ans[0] = tvd;
    ans[1] = time_sum;
    return ans;
}

void TVD(double temp, int IterMax, int R){
    int Iter[36] = {2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
                    10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000,
                    100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000,
                    1000000};

    string label[8] = {"Metropolis", "Rejection Free",
                       "Metropolis w/ AC of Size 10", "RF w/ AC of Size 10",
                       "Metropolis w/ AC of Size 50", "RF w/ AC of Size 50",
                       "Metropolis w/ AC of Size 100", "RF w/ AC of Size 100"};

    ofstream outfile;
    outfile.open("IsingAlternatingChains.txt");
    outfile << "Iter," << "TVD," << "Time," << "Algorithm" << endl;

    for (int i = 0; i < 100; i++) {
        if(Iter[i] > IterMax)
        {
            break;
        }
        cout << Iter[i] << endl;
        for (int j = 0; j < R; j++) {
            double * ans0 = Metropolis(temp, Iter[i]);
            outfile << Iter[i] << "," << ans0[0] << "," << ans0[1] << "," << label[0] << endl;

            double * ans1 = RejectionFree(temp, Iter[i]);
            outfile << Iter[i] << "," << ans1[0] << "," << ans1[1] << "," << label[1] << endl;

            double * ans2 = Metropolis_AlternatingChains(temp, Iter[i], 10);
            outfile << Iter[i] << "," << ans2[0] << "," << ans2[1] << "," << label[2] << endl;

            double * ans3 = RejectionFree_AlternatingChains(temp, Iter[i], 10);
            outfile << Iter[i] << "," << ans3[0] << "," << ans3[1] << "," << label[3] << endl;

            double * ans4 = Metropolis_AlternatingChains(temp, Iter[i], 50);
            outfile << Iter[i] << "," << ans4[0] << "," << ans4[1] << "," << label[4] << endl;

            double * ans5 = RejectionFree_AlternatingChains(temp, Iter[i], 50);
            outfile << Iter[i] << "," << ans5[0] << "," << ans5[1] << "," << label[5] << endl;

            double * ans6 = Metropolis_AlternatingChains(temp, Iter[i], 100);
            outfile << Iter[i] << "," << ans6[0] << "," << ans6[1] << "," << label[6] << endl;

            double * ans7 = RejectionFree_AlternatingChains(temp, Iter[i], 100);
            outfile << Iter[i] << "," << ans7[0] << "," << ans7[1] << "," << label[7] << endl;

        }
    }
    outfile.close();
}

int main() {
    double temp = 1;
    for(int i = 0; i < 1; i ++) {
        TVD(temp, 1000000, 100);
    }
    return 0;
}

