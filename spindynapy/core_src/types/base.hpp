#ifndef __TYPES_BASE_HPP__
#define __TYPES_BASE_HPP__

/**
 * Заголовки с базовыми интерфейсами для базовых типов/классов:
 *    - стандартные миксины
 *    - стандартные единицы абстракции (материал, спин, координаты, etc.)
 *
 * Запрещают конструктор по умолчанию.
 * Все процессы должны зависеть от этих интерфейсов.
 */

#include <format>
#include <memory>
#include <string>
#include <utility>

namespace spindynapy {

/**
 * (Абстрактный) Миксин, добавляющий методы __str__ и __repr__
 * для создаваемых объектов для удобства принтинга
 */
class StrPresentationMixin {
  public:
    StrPresentationMixin() = default;

    /**
     * Строковое (людское...) представление
     */
    virtual std::string __str__() const { return nullptr; };

    /**
     * Python-выражение объекта
     */
    virtual std::string __repr__() const { return nullptr; };

    virtual ~StrPresentationMixin() {};
};

/**
 * (Интерфейс) Базовая единица абстракции - свойства материала.
 * Хранится в регистре, в остальных - в виде ссылки, добываемой из регистра.
 */
struct IMaterial : public StrPresentationMixin {
    IMaterial() = default;
    virtual ~IMaterial() {};
};

/**
 * (Интерфейс) Базовая единица абстракции - координаты.
 * Описывает расположение в пространстве объекта.
 */
struct ICoordinates : public StrPresentationMixin {
    ICoordinates() = default;
    virtual ~ICoordinates() {};
};

/**
 * (Интерфейс) Базовая единица абстракции - направление, вектор.
 * Описывает направление сил, движения, etc.
 */
struct IDirection : public StrPresentationMixin {
    IDirection() = default;
    virtual ~IDirection() {};
};

/**
 * (Интерфейс) Структурная базовая единица абстракции - момент.
 * Магнитный, квантово-механический, etc.
 * ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.
 */
class IMoment : public StrPresentationMixin {
  public:
    IMoment() = default;
    virtual ~IMoment() {};

    /**
     * Получить вектор момента
     */
    virtual IDirection &getDirection() = 0;

    /**
     * Получить расположение момента
     */
    virtual ICoordinates &getCoordinates() = 0;
};

/**
 * (Интерфейс) Структурная единица абстракции
 * Спиновый момент. ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.
 */
class ISpin : virtual public IMoment {
  public:
    ISpin() = default;
    virtual ~ISpin() {};
};

}; // namespace spindynapy

namespace spindynapy::doc {

constexpr char module_types_base[] = "Модуль предоставляет базовые абстракные классы и интерфейсы \n"
                                     "  - базовые единицы абстракции \n"
                                     "  - миксины \n\n"
                                     "исключительно интерфейсные, не должны создаваться сами";

constexpr char StrPresentationMixin[] = "(Абстрактный) Миксин, добавляющий методы __str__ и __repr__ \n"
                                        "для создаваемых объектов для удобства принтинга";
constexpr char StrPresentationMixin___str__[] = "Строковое (людское...) представление";
constexpr char StrPresentationMixin___repr__[] = "Python-выражение объекта";

constexpr char ICoordinates[] = "(Интерфейс) Базовая единица абстракции - координаты.\n"
                                "Описывает расположение в пространстве объекта.";

constexpr char IDirection[] = "(Интерфейс) Базовая единица абстракции - направление, вектор.\n"
                              "Описывает направление сил, движения, etc.";

constexpr char IMoment[] = "(Интерфейс) Структурная базовая единица абстракции - момент.\n"
                           "Магнитный, квантово-механический, etc."
                           "ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.";
constexpr char IMoment_getDirection[] = "Получить вектор момента";
constexpr char IMoment_getCoordinates[] = "Получить расположение момента";

constexpr char ISpin[] = "(Интерфейс) Структурная единица абстракции. Спиновый момент.\n"
                         "ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.";

constexpr char IMaterial[] = "(Интерфейс) Базовая единица абстракции - свойства материала.\n"
                             "Хранится в регистре, в остальных - в виде ссылки, добываемой из регистра.\n";

}; // namespace spindynapy::doc

#endif // ! __TYPES_BASE_HPP__
