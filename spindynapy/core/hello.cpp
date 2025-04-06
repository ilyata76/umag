#include "hello.hpp"
#include <iostream>

HelloClass::HelloClass(const std::string& name) : name_(name) {}

void HelloClass::greet() const {
    std::cout << "Hello, " << name_ << "!" << std::endl;
}

/**
 * amogus
 */
void hello_function(const std::string& name) {
    std::cout << "Hello from function, " << name << "!" << std::endl;
}
