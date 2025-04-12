#ifndef __REGISTRIES_BASE_HPP__
#define __REGISTRIES_BASE_HPP__

/**
 * Регистры - контейнеры, хранящие в себе экземпляры переиспользуемых классов,
 * которые задают базовые настройки симуляции (например, свойства конкретного материала).
 * У регистров может быть расширенный функционал.
 * Реализуются для паттерна "легковес". Живут всю программу.
 */

#include "../types/base.hpp"

namespace spindynapy {

/**
 * Базовый интерфейс регистра-контейнера
 */
class IRegistry : public StrPresentationMixin {
  public:
    IRegistry() = default;
    virtual ~IRegistry() = 0;
};

/**
 * Интерфейс для регистра-контейнера для материалов
 */
class IMaterialRegistry : public IRegistry {
  public:
    IMaterialRegistry() = default;
    virtual ~IMaterialRegistry() = 0;
};

/**
 * Интерфейс для регистра-контейнера для разных областей (регионов геометрии) образца
 */
class IRegionRegistry : public IRegistry {
  public:
    IRegionRegistry() = default;
    virtual ~IRegionRegistry() = 0;
};

}; // namespace spindynapy

namespace spindynapy::doc {

constexpr char module_registries[] =
    " Регистры - контейнеры, хранящие в себе экземпляры переиспользуемых классов,\n"
    " которые задают базовые настройки симуляции (например, свойства конкретного материала).\n"
    " У регистров может быть расширенный функционал.\n"
    " Реализуются для паттерна \"легковес\". Живут всю программу.";

constexpr const char *module_registries_base = module_registries;

constexpr char IRegistry[] = "Базовый интерфейс регистра-контейнера";
constexpr char IMaterialRegistry[] = "Интерфейс для регистра-контейнера для материалов";
constexpr char IRegionRegistry[] = "Интерфейс для регистра-контейнера для разных областей материалов";

}; // namespace spindynapy::doc

#endif // ! __REGISTRIES_BASE_HPP__
