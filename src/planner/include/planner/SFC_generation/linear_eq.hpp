#ifndef LINE_EQUATIONS_HPP
#define LINE_EQUATIONS_HPP

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <Eigen/Dense>

using namespace std;

namespace line_eq
{

class LineEquations
{
public:
    using Point = std::pair<int, int>;

    LineEquations() = default;

    LineEquations(const std::vector<Point>& points,
                  const Point& seed)
    {   
        //cout<<"seed = "<<seed.first<<" "<<seed.second<<" "<<endl;
        //cout<<"points = "<<points<<endl;
        setPoints(points);
        setSeed(seed);
    }

    // -----------------------------
    // Setters
    // -----------------------------
    void setPoints(const std::vector<Point>& points)
    {
        points_ = points;
        computed_ = false;
    }

    void setSeed(const Point& seed)
    {
        seed_ = seed;
        has_seed_ = true;
        computed_ = false;
    }

    // ----------------------------------------------------
    // Compute Ax <= b from polygon edges
    // Inequality is oriented so SEED satisfies it
    // ----------------------------------------------------
    void compute()
    {
        equations_.clear();

        const int n = static_cast<int>(points_.size());
        if (n < 2 || !has_seed_)
        {
            A_.resize(0, 2);
            b_.resize(0);
            computed_ = true;
            return;
        }



        double area = 0.0;
        for (int i = 0; i < n; ++i) {
            auto& p = points_[i];
            auto& q = points_[(i + 1) % n];
            area += (p.first * q.second - q.first * p.second);
        }
        if (area < 0.0) {
            std::reverse(points_.begin(), points_.end());
        }

        A_.resize(n, 2);
        b_.resize(n);

        const double sx = seed_.first;
        const double sy = seed_.second;

        // ------------------------------------------------------------------
        // 2. Build half-spaces
        // ------------------------------------------------------------------

        for (int i = 0; i < n; ++i)
        {
            double x1 = points_[i].first;
            double y1 = points_[i].second;

            double x2 = points_[(i + 1) % n].first;
            double y2 = points_[(i + 1) % n].second;

            // Standard line coefficients
            double a = y1 - y2;
            double b = x2 - x1;
            double c = x1 * y2 - x2 * y1;

            // Evaluate line at seed
            double eval = a * sx + b * sy + c;

            //cout<<sx<<" "<<sy<<" "<<a<<" "<<b<<" "<<c<<eval<<endl;



            // ❗ Flip inequality if seed is on wrong side
            constexpr double EPS = 1e-6;
            if (eval > EPS)
            {

                //cout<<"yes_its more"<<a<<" "<<b<<" "<<c<<endl;
                a = -a;
                b = -b;
                c = -c;
            }

            //cout<<a<<" "<<b<<" "<<c<<endl;

            // Store Ax <= b
            A_(i, 0) = a;
            A_(i, 1) = b;
            b_(i)    = -c;

            // Human-readable inequality
            equations_.emplace_back(
                std::to_string(a) + "x + " +
                std::to_string(b) + "y <= " +
                std::to_string(-c)
            );
        }

            // ------------------------------------------------------------------
            // 3. Final safety check: seed must satisfy all constraints
            // ------------------------------------------------------------------
            constexpr double EPS = 1e-6;

            for (int i = 0; i < n; ++i) {
                double lhs = A_(i,0) * sx + A_(i,1) * sy;
                if (lhs > b_(i) + EPS) {
                    std::cerr << "[ERROR] Seed violates constraint " << i << "\n";
                    std::cerr << "  lhs = " << lhs
                            << ", rhs = " << b_(i) << std::endl;
                    std::abort();   // or throw
                }
            }

        computed_ = true;
    }

    // -----------------------------
    // Accessors
    // -----------------------------
    const Eigen::MatrixXd& A() const { return A_; }
    const Eigen::VectorXd& b() const { return b_; }
    const std::vector<std::string>& equations() const { return equations_; }

    // -----------------------------
    // Debug print
    // -----------------------------
    void print(std::ostream& os = std::cout) const
    {
        if (!computed_)
            os << "Warning: compute() not called yet.\n";

        os << "Seed voxel : (" << seed_.first << ", " << seed_.second << ")\n\n";

        os << "Matrix A:\n" << A_ << "\n\n";
        os << "Vector b:\n" << b_ << "\n\n";

        os << "Inequalities (Ax <= b):\n";
        for (const auto& eq : equations_)
            os << eq << "\n";
    }

private:
    std::vector<Point> points_;
    Point seed_{0, 0};
    bool has_seed_ = false;

    Eigen::MatrixXd A_;
    Eigen::VectorXd b_;
    std::vector<std::string> equations_;

    bool computed_ = false;
};

} // namespace line_eq

#endif // LINE_EQUATIONS_HPP
