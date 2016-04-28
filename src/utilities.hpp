#pragma once

#include <iostream>
#include <sstream>
#include <random>

namespace hmc { namespace utilities {
    std::default_random_engine random_generator;

    void set_random_state_true_random()
    {
        random_generator.seed(std::random_device{}());
    }

    std::string get_random_state()
    {
        std::ostringstream s;
        s << random_generator;
        return s.str();
    }

    void set_random_state(std::string state)
    {
        std::istringstream {state} >> random_generator;
    }

    double uniform(double a, double b)
    {
        return std::uniform_real_distribution<double>(a, b)(random_generator);
    }
}}
