#ifndef HELLO_HPP
#define HELLO_HPP

#include <string>

class HelloClass {
public:
    HelloClass(const std::string& name);
    void greet() const;
private:
    std::string name_;
};

void hello_function(const std::string& name);

#endif // !HELLO_HPP
