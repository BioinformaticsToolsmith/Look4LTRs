#include <memory>
#include <iostream>

int main() {
    auto a = std::make_shared<int>(5);
    std::cout << &a << std::endl;


}