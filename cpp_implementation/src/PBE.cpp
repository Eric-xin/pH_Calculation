#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

class PBE_Acid {
   public:
    PBE_Acid(const std::vector<double>& Ka, const std::vector<double>& pKa, int proton, int proton_ref, double conc)
        : proton_ref(proton_ref), conc(conc) {
        if (Ka.empty() && pKa.empty()) {
            throw std::invalid_argument("You must define either Ka or pKa values.");
        }
        if (proton == 0) {
            throw std::invalid_argument("The maximum proton for this acid must be defined.");
        }
        // if (proton_ref == 0) {
        //     throw std::invalid_argument("The reference proton for PBE must be defined.");
        // }

        if (Ka.empty()) {
            this->pKa = pKa;
            std::sort(this->pKa.begin(), this->pKa.end());
            this->Ka.resize(this->pKa.size());
            std::transform(this->pKa.begin(), this->pKa.end(), this->Ka.begin(), [](double pKa) { return std::pow(10, -pKa); });
        } else {
            this->Ka = Ka;
            std::sort(this->Ka.begin(), this->Ka.end(), std::greater<double>());
            this->pKa.resize(this->Ka.size());
            std::transform(this->Ka.begin(), this->Ka.end(), this->pKa.begin(), [](double Ka) { return -std::log10(Ka); });
        }

        _Ka_temp.push_back(1.0);
        _Ka_temp.insert(_Ka_temp.end(), this->Ka.begin(), this->Ka.end());

        for (int i = 0; i <= this->Ka.size(); ++i) {
            this->proton_vector.push_back(proton - i - proton_ref);
        }
    }

    std::vector<double> alpha(double pH) const {
        double h3o = std::pow(10, -pH);
        std::vector<double> h3o_pow(_Ka_temp.size());
        std::vector<double> Ka_prod(_Ka_temp.size());

        for (size_t i = 0; i < _Ka_temp.size(); ++i) {
            h3o_pow[i] = std::pow(h3o, _Ka_temp.size() - 1 - i);
        }

        std::partial_sum(_Ka_temp.begin(), _Ka_temp.end(), Ka_prod.begin(), std::multiplies<double>());

        std::vector<double> h3o_Ka(_Ka_temp.size());
        for (size_t i = 0; i < _Ka_temp.size(); ++i) {
            h3o_Ka[i] = h3o_pow[i] * Ka_prod[i];
        }

        double den = std::accumulate(h3o_Ka.begin(), h3o_Ka.end(), 0.0);
        std::vector<double> result(_Ka_temp.size());
        std::transform(h3o_Ka.begin(), h3o_Ka.end(), result.begin(), [den](double val) { return val / den; });

        return result;
    }

    int get_proton(int index) const {
        return proton_vector[index];
    }

    double get_conc() const {
        return conc;
    }

    size_t get_proton_vector_size() const {
        return proton_vector.size();
    }

   private:
    std::vector<double> Ka;
    std::vector<double> pKa;
    std::vector<double> _Ka_temp;
    std::vector<int> proton_vector;
    int proton_ref;
    double conc;
};

class PBE_calc {
   public:
    PBE_calc(const std::vector<PBE_Acid>& acids, double Kw = 1.01e-14)
        : acids(acids), Kw(Kw) {}

    double PBE_error(double pH) const {
        double h3o = std::pow(10, -pH);
        double oh = Kw / h3o;
        double P_error = h3o - oh;

        for (const auto& acid : acids) {
            auto alpha = acid.alpha(pH);
            for (size_t i = 0; i < alpha.size(); ++i) {
                P_error += acid.get_conc() * acid.get_proton(i) * alpha[i];
            }
        }

        return std::abs(P_error);
    }

    double pH_calc(double guess = 7.0, bool guess_est = false, int est_num = 1500, double tol = 1e-5) {
        if (guess_est) {
            std::vector<double> phs(est_num);
            double step = 14.0 / (est_num - 1);
            for (int i = 0; i < est_num; ++i) {
                phs[i] = i * step;
            }

            double min_diff = std::numeric_limits<double>::max();
            for (double ph : phs) {
                double diff = PBE_error(ph);
                if (diff < min_diff) {
                    min_diff = diff;
                    guess = ph;
                }
            }
        }

        double pH = guess;
        double diff = PBE_error(pH);
        double learning_rate = 0.001;  // Stochastic Gradient Descent learning rate
        int max_iterations = 10000;    // Maximum number of iterations to prevent infinite loop

        for (int i = 0; i < max_iterations && diff > tol; ++i) {
            double gradient = (PBE_error(pH + tol) - diff) / tol;  // Numerical gradient
            pH -= learning_rate * gradient;                        // Update pH using gradient descent
            diff = PBE_error(pH);                                  // Recalculate the difference
        }

        if (diff > tol) {
            throw std::runtime_error("Failed to converge to the desired tolerance.");
        }

        return pH;
    }

   private:
    std::vector<PBE_Acid> acids;
    double Kw;
};

// int main() {
//     // std::vector<double> Ka = {1.0e-3, 1.0e-5, 1.0e-7};
//     // std::vector<double> pKa = {3.0, 5.0, 7.0};
//     // Acid acetic(Ka, {}, 1, 0.1);
//     // Acid formic({}, pKa, 1, 0.1);
//     // std::vector<Acid> species = {acetic, formic};

//     std::vector<double> pKa_p = {1.97, 6.82, 12.5};
//     PBE_Acid phosphoric({}, pKa_p, 3, 1, 0.01);
//     PBE_Acid NH4({}, {9.25}, 1, 1, 0.01*2);
//     std::vector<PBE_Acid> species = {phosphoric, NH4};
//     // std::vector<Acid> species = {NH4};

//     PBE_calc cbe(species);
//     double pH = cbe.pH_calc(7.0, true, 1500, 1e-5);
//     std::cout << "The pH is: " << pH << std::endl;

//     return 0;
// }

int main() {
    printf("The program is running...\n");
    printf("Program name: Calculation of pH with Proton Balance Equation (PBE)\n");
    printf("Author: Eric Xin\n");
    printf("----------------------------------------\n");

    printf("Please enter the number of components: ");
    int num_components;
    scanf("%d", &num_components);

    printf("Confirmed, the number of components is %d.\n", num_components);
    printf("----------------------------------------\n");
    printf("Now please enter the information for each component.\n");
    std::vector<PBE_Acid> acids;
    while (num_components) {
        printf("Using Ka or pKa values? (1 for Ka, 2 for pKa): ");
        int sub_choice;
        scanf("%d", &sub_choice);

        if (sub_choice == 1) {
            printf("Please enter the number of Ka values: ");
            int num_Ka;
            scanf("%d", &num_Ka);
            std::vector<double> Ka(num_Ka);
            printf("Please enter the Ka values: ");
            for (int i = 0; i < num_Ka; ++i) {
                scanf("%lf", &Ka[i]);
            }

            printf("Please enter the maximum proton of the acid: ");
            int proton;
            scanf("%d", &proton);

            printf("Please enter the reference proton for PBE: ");
            int proton_ref;
            scanf("%d", &proton_ref);

            printf("Please enter the concentration of the acid: ");
            double conc;
            scanf("%lf", &conc);

            acids.push_back(PBE_Acid(Ka, {}, proton, proton_ref, conc));
        } else if (sub_choice == 2) {
            printf("Please enter the number of pKa values: ");
            int num_pKa;
            scanf("%d", &num_pKa);
            std::vector<double> pKa(num_pKa);
            printf("Please enter the pKa values: ");
            for (int i = 0; i < num_pKa; ++i) {
                scanf("%lf", &pKa[i]);
            }

            printf("Please enter the maximum proton of the acid: ");
            int proton;
            scanf("%d", &proton);

            printf("Please enter the reference proton for PBE: ");
            int proton_ref;
            scanf("%d", &proton_ref);

            printf("Please enter the concentration of the acid: ");
            double conc;
            scanf("%lf", &conc);

            acids.push_back(PBE_Acid({}, pKa, proton, proton_ref, conc));
        }

        num_components--;
    }

    printf("----------------------------------------\n");
    printf("Please enter the value of Kw: (input 0 for default value 1.01e-14): ");
    double Kw;
    scanf("%lf", &Kw);
    if (Kw == 0) {
        Kw = 1.01e-14;
        printf("Using default value of Kw: 1.01e-14\n");
    } else {
        printf("The value of Kw is: %.2e\n", Kw);
    }
    printf("----------------------------------------\n");
    int est = 1;

    PBE_calc pbe(acids, Kw);
    double pH;
    if (est) {
        printf("All parameters are set. Estimating the initial pH...\n");
        pH = pbe.pH_calc(7.0, true, 1500, 1e-5);
    } else {
        printf("Please enter the initial guess of pH: (or input 0 for default value 7.0): ");
        double guess;
        scanf("%lf", &guess);
        if (guess == 0) {
            guess = 7.0;
            printf("Using default value of pH: 7.0\n");
        } else {
            printf("The initial guess of pH is: %.2f\n", guess);
        }
        printf("All parameters are set. Calculating the pH...\n");
        pH = pbe.pH_calc(guess, false, 1500, 1e-5);
    }
    printf("----------------------------------------\n");
    printf("The pH is: %.5f\n", pH);
    printf("----------------------------------------\n");
    printf("Do you want a calculation report? (1 for yes, 0 for no): ");
    int report;
    scanf("%d", &report);
    if (report) {
        printf("Here is a digest of all the components:\n");
        printf("----------------------------------------\n");
        for (const auto& acid : acids) {
            printf("For component No. %d:\n", (int)(&acid - &acids[0]) + 1);
            printf("Component with concentration %.2e\n", acid.get_conc());
            printf("Alpha values: at equilibrium pH = %.5f\n", pH);
            auto alpha = acid.alpha(pH);
            for (size_t i = 0; i < alpha.size(); ++i) {
                printf("Alpha_%d: %.5f\n", acid.get_proton(i), alpha[i]);
            }
            printf("----------------------------------------\n");
        }
    }
    printf("Thank you for using the program. Goodbye!\n");
}