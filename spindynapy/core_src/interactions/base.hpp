#ifndef __INTERACTIONS_BASE_HPP__
#define __INTERACTIONS_BASE_HPP__

/**
 * Интерфейсы взаимодействий между элементами системы.
 * Таковыми могут быть потенциальные поля, силы, etc.
 */

#include "../types/base.hpp"

namespace spindynapy {

/**
 * Базовый интерфейс взаимодействий (сил, полей, etc.)
 */
class IInteraction : public StrPresentationMixin {
  public:
    IInteraction() = default;
    virtual ~IInteraction() = 0;
};

}; // namespace spindynapy

namespace spindynapy::doc {

constexpr char module_interactions[] = "Интерфейсы взаимодействий между элементами системы.\n"
                                       "Таковыми могут быть потенциальные поля, силы, etc. \n";

constexpr const char *module_interactions_base = module_interactions;

constexpr char IInteraction[] = "Базовый интерфейс взаимодействий (сил, полей, etc.)";

}; // namespace spindynapy::doc

#endif // ! __INTERACTIONS_BASE_HPP__
