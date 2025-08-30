#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

// Class definitions
class Vessel {
    public:
        Vessel(double m, double A, double K, double y0, double v0, double h0)
            : m(m), A(A), K(K), y0(y0), v0(v0), h0(h0) {}; // Constructor
        ~Vessel() {
            // Nothing to do for this class
        }; // Destructor

        // Getters
        double getMass();
        double getArea();
        double getFlowRateCoefficient();
        double getInitialHeight();
        double getInitialVelocity();
        double getInitialWaterLevel();

        // Setters
        void setMass(double m);
        void setArea(double A);
        void setFlowRateCoefficient(double K);
        void setInitialHeight(double y0);
        void setInitialVelocity(double v0);
        void setInitialWaterLevel(double h0);

        // Numerical methods
        vector<vector<double>> eulerMethod(double t, double dt, double g, double rho, double Cd);
        vector<vector<double>> rungeKuttaMethod(double t, double dt, double g, double rho, double Cd);

    private:
        double m, A, K, y0, v0, h0;
};

double fv(double m, double A, double K, double y_prev, double v_prev, double h_prev, double g, double rho, double Cd) {
    return ((rho * A) / (m + rho * A * h_prev)) * (2 * g * h_prev - g * y_prev - 0.5 * Cd * v_prev * abs(v_prev));
};

double fh(double m, double A, double K, double y_prev, double h_prev, double rho) {
    return ((K / (rho * A)) * ((m / (rho * A)) + y_prev - h_prev));
};


vector<vector<double>> Vessel::eulerMethod(double t, double dt, double g, double rho, double Cd) {
    const int length = t / dt + 1;
    vector<vector<double>> results(3, vector<double>(length));
    // Initial conditions
    results[0][0] = y0;
    results[1][0] = v0;
    results[2][0] = h0;

    for (int i = 1; i < length; i++)
    {
        // Temporary variables
        double y_prev = results[0][i - 1]; // While these temporary variables are not necessary,
        double v_prev = results[1][i - 1]; // they make the code more readable
        double h_prev = results[2][i - 1];

        // Solving for y(t)
        results[0][i] = y_prev + dt * v_prev;

        // Solving for v(t)
        results[1][i] = v_prev + dt * fv(m, A, K, y_prev, v_prev, h_prev, g, rho, Cd);

        // Solving for h(t)
        results[2][i] = h_prev + dt * fh(m, A, K, y_prev, h_prev, rho);
    }

    return results;
}

vector<vector<double>> Vessel::rungeKuttaMethod(double t, double dt, double g, double rho, double Cd) {
    const int length = t / dt + 1;
    vector<vector<double>> results(3, vector<double>(length));

    // Initial conditions
    results[0][0] = y0;
    results[1][0] = v0;
    results[2][0] = h0;

    for (int i = 1; i < length; i++) {
        // Previous values
        double y_prev = results[0][i - 1];
        double v_prev = results[1][i - 1];
        double h_prev = results[2][i - 1];

        // k1 values
        double k1_y = v_prev;
        double k1_v = fv(m, A, K, y_prev, v_prev, h_prev, g, rho, Cd);
        double k1_h = fh(m, A, K, y_prev, h_prev, rho);

        // k2 values
        double k2_y = v_prev + 0.5 * dt * k1_v;
        double k2_v = fv(m, A, K, y_prev + 0.5 * dt * k1_y, v_prev + 0.5 * dt * k1_v, h_prev + 0.5 * dt * k1_h, g, rho, Cd);
        double k2_h = fh(m, A, K, y_prev + 0.5 * dt * k1_y, h_prev + 0.5 * dt * k1_h, rho);

        // k3 values
        double k3_y = v_prev + 0.5 * dt * k2_v;
        double k3_v = fv(m, A, K, y_prev + 0.5 * dt * k2_y, v_prev + 0.5 * dt * k2_v, h_prev + 0.5 * dt * k2_h, g, rho, Cd);
        double k3_h = fh(m, A, K, y_prev + 0.5 * dt * k2_y, h_prev + 0.5 * dt * k2_h, rho);

        // k4 values
        double k4_y = v_prev + dt * k3_v;
        double k4_v = fv(m, A, K, y_prev + dt * k3_y, v_prev + dt * k3_v, h_prev + dt * k3_h, g, rho, Cd);
        double k4_h = fh(m, A, K, y_prev + dt * k3_y, h_prev + dt * k3_h, rho);

        // Update y(t), v(t), and h(t)
        results[0][i] = y_prev + (dt / 6.0) * (k1_y + 2 * k2_y + 2 * k3_y + k4_y);
        results[1][i] = v_prev + (dt / 6.0) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
        results[2][i] = h_prev + (dt / 6.0) * (k1_h + 2 * k2_h + 2 * k3_h + k4_h);
    }

    return results;
}

int main() {
    string line;
    ifstream inpFile("parameters.txt", ios::in);

    int scheme;
    double T, dt, g, rho, Cd;

    if (inpFile.good()) {
        // Initialize vector of vessels
        vector<Vessel> vessels;


        while (getline(inpFile, line)) {
            // Skip lines that start with '#'
            if (line[0] == '#') {
                continue;
            }

            stringstream ss(line);
            if (!(ss >> scheme >> T >> dt >> g >> rho >> Cd)) {
                cerr << "Error: Failed to parse the second line of parameters.txt. "
                    << "Ensure it contains valid numeric values." << endl;
                cerr << "Line content: " << line << endl;
                return 1;
            }

            while (getline(inpFile, line)) {
                if (line[0] == '#') continue;
                stringstream ss(line);
                double m, A, K, y0, v0, h0;
                ss >> m >> A >> K >> y0 >> v0 >> h0;
                vessels.push_back(Vessel(m, A, K, y0, v0, h0));
            }
        }   

        ofstream outFile("output.txt", ios::out | ios::trunc);
        if (!outFile.is_open()) {
            cerr << "Error: Could not open output.txt for writing!" << endl;
            return 1;
        }

        int numTimeSteps = T/dt + 1;

        vector<vector<vector<double>>> allResults;
        for (auto& vessel : vessels) {
            vector<vector<double>> results = (scheme == 0)
                ? vessel.eulerMethod(T, dt, g, rho, Cd)
                : vessel.rungeKuttaMethod(T, dt, g, rho, Cd);

            allResults.push_back(results);
        }

        for (int t = 0; t < numTimeSteps; t++) {
            double currentTime = t * dt;

            outFile << setw(15) << currentTime;

            for (const auto& results : allResults) {

                outFile << setw(15) << results[0][t]
                        << setw(15) << results[1][t]
                        << setw(15) << results[2][t];
            }
            outFile << endl;
        }

        cout << "Output written to output.txt" << endl;
    } else {
        cerr << "Error: Could not open parameters.txt!" << endl;
    }

    return 0;
}
