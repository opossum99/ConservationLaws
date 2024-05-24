#ifndef CONSERVATIONLAWSPROJECT_GEOMETRY_HPP
#define CONSERVATIONLAWSPROJECT_GEOMETRY_HPP

#include <vector>
#include <functional>
#include <set>
#include <map>

class Obstacle {
public:
    Obstacle(double x0, double y0, double length, double height);

private:
    double x0_;
    double y0_;
    double length_;
    double height_;
};

class Geometry {
public:
    explicit Geometry(int N_x, int N_y);

    explicit Geometry(std::vector<Obstacle>);

    const std::set<int> &interior() { return indexes_inner_; }

    const std::map<std::pair<std::string, std::string>, std::set<int>> &ghost() { return indexes_ghost_; }

    const std::set<int>& ghost(const std::pair<std::string, std::string>& str) { return indexes_ghost_[str]; }

private:
    std::vector<std::pair<double, double>> xy_;
    std::set<int> indexes_inner_;
    std::map<std::pair<std::string, std::string>, std::set<int>> indexes_ghost_;
};


#endif //CONSERVATIONLAWSPROJECT_GEOMETRY_HPP
