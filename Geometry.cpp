#include "Geometry.hpp"

Geometry::Geometry(int N_x, int N_y) {
    for (int i = 0; i < N_x * N_y; i++) {
        if (i < N_x) {
            indexes_ghost_[{"down", "firm"}].insert(i);
        } else if (i >= N_x * (N_y - 1)){
            indexes_ghost_[{"up", "free"}].insert(i);
        } else if (i % N_x == 0){
            indexes_ghost_[{"left", "firm"}].insert(i);
        } else if (i % N_x == N_x - 1){
            indexes_ghost_[{"right", "free"}].insert(i);
        } else {
            indexes_inner_.insert(i);
        }
    }
}
