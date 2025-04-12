#ifndef __SOLVERS_BASE_HPP__
#define __SOLVERS_BASE_HPP__

/**
 * Решатели (интеграторы) берут на себя задачу эволюционирования
 * или любого другого прогрессирования системы в зависимости от
 * предоставленных внешних параметров.
 */

#include "../types/base.hpp"

namespace spindynapy {

/**
 * Базовый интерфейс решателя-интегратора
 */
class ISolver {
  public:
    ISolver() = default;
    virtual ~ISolver() = 0;

    virtual std::string __str__() const { return nullptr; };
    virtual std::string __repr__() const { return nullptr; };
};

}; // namespace spindynapy

namespace spindynapy::doc {

constexpr char module_solvers[] = "Модуль решателей\n"
                                  "Решатели (интеграторы) берут на себя задачу эволюционирования\n"
                                  "или любого другого прогрессирования системы в зависимости от\n"
                                  "предоставленных внешних параметров.";
constexpr const char *module_solvers_base = module_solvers;

constexpr char ISolver[] = "Базовый интерфейс решателя-интегратора";

}; // namespace spindynapy::doc

#endif // ! __SOLVERS_BASE_HPP__
