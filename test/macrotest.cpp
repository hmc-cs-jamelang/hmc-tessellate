#define PRINTEM(...) __VA_ARGS__
#include <iostream>

int main() {
    std::cout << __LINE__ << std::endl;
    PRINTEM(
        std::cout << __LINE__ << std::endl;
        std::cout << __LINE__ << std::endl;
        std::cout << __LINE__ << std::endl;
        std::cout << __LINE__ << std::endl;
        std::cout << __LINE__ << std::endl;
        std::cout << __LINE__ << std::endl;
    )
    std::cout << __LINE__ << std::endl;
    return 0;
}
