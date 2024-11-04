// #include <iostream>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

class Acid {
   public:
    Acid(const std::vector<double> &Ka, const std::vector<double> &pKa, int charge, double conc)
        : charge(charge), conc(conc) {
        if (Ka.empty() && pKa.empty()) {
            throw std::invalid_argument("You must define either Ka or pKa values.");
        }
        // if (charge == 0) {
        //     throw std::invalid_argument("The maximum charge for this acid must be defined.");
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
            this->charge_vector.push_back(charge - i);
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

    inline void print_acid_data() const {
        if (pKa.empty()) {
            printf("Ka values: ");
            for (const auto &val : Ka) {
                printf("%.2e ", val);
            }
            printf("\n");
        } else {
            printf("pKa values: ");
            for (const auto &val : pKa) {
                printf("%.2f ", val);
            }
            printf("\n");
        }

        printf("Charge: %d\n", charge);
        printf("Concentration: %.2e\n", conc);
    }

    double get_conc() const {
        return conc;
    }

    const std::vector<int> &get_charge_vector() const {
        return charge_vector;
    }

   private:
    std::vector<double> Ka;
    std::vector<double> pKa;
    std::vector<double> _Ka_temp;
    std::vector<int> charge_vector;
    int charge;
    double conc;
};

class CBE_calc {
   public:
    CBE_calc(const std::vector<Acid> &species, double Kw = 1.01e-14)
        : species(species), Kw(Kw) {}

    double Charge_diff(double pH) const {
        double h3o = std::pow(10, -pH);
        double oh = Kw / h3o;
        double x = h3o - oh;

        for (const auto &s : species) {
            auto alpha = s.alpha(pH);
            for (size_t i = 0; i < alpha.size(); ++i) {
                x += s.get_conc() * s.get_charge_vector()[i] * alpha[i];
            }
        }

        return std::abs(x);
    }

    double pH_calc(double guess = 7.0, bool guess_est = false, int est_num = 1500, double tol = 1e-5) {
        printf("Calculating pH...\n");
        if (guess_est) {
            std::vector<double> phs(est_num);
            double step = 14.0 / (est_num - 1);
            for (int i = 0; i < est_num; ++i) {
                phs[i] = i * step;
            }

            double min_diff = std::numeric_limits<double>::max();
            for (double ph : phs) {
                double diff = Charge_diff(ph);
                if (diff < min_diff) {
                    min_diff = diff;
                    guess = ph;
                }
            }
        }

        double pH = guess;
        double diff = Charge_diff(pH);
        double learning_rate = 0.001;  // Stochastic Gradient Descent learning rate
        int max_iterations = 10000;    // Maximum number of iterations to prevent infinite loop

        for (int i = 0; i < max_iterations && diff > tol; ++i) {
            double gradient = (Charge_diff(pH + tol) - diff) / tol;  // Numerical gradient
            pH -= learning_rate * gradient;                          // Update pH using gradient descent
            diff = Charge_diff(pH);                                  // Recalculate the difference

            // printf("Iteration: %d, pH: %.5f, diff: %.5e, gradient: %.5e\n", i, pH, diff, gradient);
        }

        if (diff > tol) {
            throw std::runtime_error("Failed to converge to the desired tolerance.");
        }

        return pH;
    }

   private:
    std::vector<Acid> species;
    double Kw;
};

// int main() {
//     // std::vector<double> Ka = {1.0e-3, 1.0e-5, 1.0e-7};
//     // std::vector<double> pKa = {3.0, 5.0, 7.0};
//     // Acid acetic(Ka, {}, 1, 0.1);
//     // Acid formic({}, pKa, 1, 0.1);
//     // std::vector<Acid> species = {acetic, formic};

//     std::vector<double> pKa_p = {1.97, 6.82, 12.5};
//     Acid phosphoric({}, pKa_p, 0, 0.01);
//     Acid NH4({}, {9.25}, 1, 0.01*3);
//     std::vector<Acid> species = {phosphoric, NH4};
//     // std::vector<Acid> species = {NH4};

//     CBE_calc cbe(species);
//     double pH = cbe.pH_calc(7.0, true, 1500, 1e-5);
//     std::cout << "The pH is: " << pH << std::endl;

//     return 0;
// }

int main() {
    printf("The program is running...\n");
    printf("Program name: Calculation of pH with Charge Balance Equation (CBE)\n");
    printf("Author: Eric Xin\n");
    printf("----------------------------------------\n");

    printf("Please enter the number of components: ");
    int num_components;
    scanf("%d", &num_components);

    printf("Confirmed, the number of components is %d.\n", num_components);
    printf("----------------------------------------\n");
    printf("Now please enter the information for each component.\n");
    std::vector<Acid> s;
    while (num_components) {
        printf("Using Ka or pKa values? (1 for Ka, 2 for pKa): ");
        int choice;
        scanf("%d", &choice);

        if (choice == 1) {
            printf("Please enter the number of Ka values: ");
            int num_Ka;
            scanf("%d", &num_Ka);
            std::vector<double> Ka(num_Ka);
            printf("Please enter the Ka values: ");
            for (int i = 0; i < num_Ka; ++i) {
                scanf("%lf", &Ka[i]);
            }

            printf("Please enter the charge of the acid: ");
            int charge;
            scanf("%d", &charge);

            printf("Please enter the concentration of the acid: ");
            double conc;
            scanf("%lf", &conc);

            s.push_back(Acid(Ka, {}, charge, conc));
        } else if (choice == 2) {
            printf("Please enter the number of pKa values: ");
            int num_pKa;
            scanf("%d", &num_pKa);
            std::vector<double> pKa(num_pKa);
            printf("Please enter the pKa values: ");
            for (int i = 0; i < num_pKa; ++i) {
                scanf("%lf", &pKa[i]);
            }

            printf("Please enter the charge of the acid: ");
            int charge;
            scanf("%d", &charge);

            printf("Please enter the concentration of the acid: ");
            double conc;
            scanf("%lf", &conc);

            s.push_back(Acid({}, pKa, charge, conc));
        } else {
            printf("Invalid choice. Please try again.\n");
            continue;
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
    // printf("Do you want to estimate the initial pH? (1 for yes, 0 for no): ");
    // int est;
    // scanf("%d", &est);
    int est = 1;

    CBE_calc cbe(s, Kw);
    double pH;
    if (est) {
        printf("All parameters are set. Estimating the initial pH...\n");
        pH = cbe.pH_calc(7.0, true, 1500, 1e-5);
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
        pH = cbe.pH_calc(guess, false, 1500, 1e-5);
    }
    printf("----------------------------------------\n");
    printf("The pH is: %.5f\n", pH);
    printf("----------------------------------------\n");
    // printf("Here is a digest of all the components:\n");
    // for (const auto& acid : s) {
    //     printf("Component with charge %d and concentration %.2e\n", acid.get_charge_vector()[0], acid.get_conc());
    //     // print the alpha values
    //     printf("Alpha values: at equilibrium pH = %.5f\n", pH);
    //     auto alpha = acid.alpha(pH);
    //     for (size_t i = 0; i < alpha.size(); ++i) {
    //         printf("Alpha_%d: %.5f\n", acid.get_charge_vector()[i], alpha[i]);
    //     }
    // }
    // printf("----------------------------------------\n");
    printf("Do you want a calculation report? (1 for yes, 0 for no): ");
    int report;
    scanf("%d", &report);
    if (report) {
        printf("Here is a digest of all the components:\n");
        printf("----------------------------------------\n");
        for (const auto &acid : s) {
            printf("For component No. %d:\n", (int)(&acid - &s[0]) + 1);
            acid.print_acid_data();
            printf("Component with charge %d and concentration %.2e\n", acid.get_charge_vector()[0], acid.get_conc());
            // print the alpha values
            printf("Alpha values: at equilibrium pH = %.5f\n", pH);
            auto alpha = acid.alpha(pH);
            for (size_t i = 0; i < alpha.size(); ++i) {
                printf("Alpha_%d: %.5f\n", acid.get_charge_vector()[i], alpha[i]);
            }
            printf("----------------------------------------\n");
        }
        // printf("----------------------------------------\n");
    }
    printf("Thank you for using the program. Goodbye!\n");
}