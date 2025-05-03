#ifndef __REGISTRIES_HPP__
#define __REGISTRIES_HPP__

/**
 * Регистры - потокобезопасные контейнеры, хранящие в себе экземпляры переиспользуемых классов,
 *   которые задают базовые настройки симуляции (например, свойства конкретного материала).
 *
 * У регистров может быть расширенный функционал.
 * "Легковесы". Живут всю программу.
 *
 * Кэши - потокобезопасные контейнеры, которые по реализации похожи на регистры,
 *    но предназначены для хранения промежуточных данных, кэш.
 */

#include "constants.hpp"
#include "types.hpp"

#include <memory>
#include <mutex>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <shared_mutex>
#include <stdexcept>
#include <string>

namespace py = pybind11;

namespace std {

/// Комбинирует существующее значение хеша `seed` с новым ключом `key`.
/// Используется для построения хешей составных объектов, таких как пары или структуры.
/// Алгоритм основан на реализации Boost::hash_combine и снижает вероятность коллизий.
template <typename T> void hash_combine(size_t &seed, T const &key) {
    hash<T> hasher;
    seed ^= hasher(key) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/// Специализация std::hash для std::pair<T1, T2>.
/// Выполняет несимметричное хеширование пары: порядок элементов имеет значение.
/// Это значит, что hash({a, b}) != hash({b, a}), если a != b.
template <typename T1, typename T2> struct hash<pair<T1, T2>> {
    std::size_t operator()(const std::pair<T1, T2> &p) const {
        std::size_t seed = 0;
        hash_combine(seed, p.first);
        hash_combine(seed, p.second);
        return seed;
    }
};

} // namespace std

namespace PYTHON_API spindynapy {

// тип для внутрянки регистра: хранимые элементы по регистровым номерам
template <typename Element>
using RegistryContainer = PYTHON_API std::unordered_map<regnum, std::shared_ptr<Element>>;

/**
 * Регистр-контейнер
 *   Предоставляет доступ к элементам, хранящимся в регистре.
 *  БЕЗ ЗАПИСИ => ПОТОКОБЕЗОПАСНЫЙ
 */
template <typename Element> class PYTHON_API Registry {
  protected:
    // внутренний контейнер для хранения элементов
    RegistryContainer<Element> _container;

  public:
    // пустой конструктор
    PYTHON_API Registry() : _container() {};
    // конструктор с предоставленным контейнером
    PYTHON_API Registry(RegistryContainer<Element> container) : _container(container) {}
    // деструктор
    ~Registry() = default;

    // строковое представление для принтинга (TODO заменить на ostream/stringstream)
    PYTHON_API std::string __str__() const {
        std::stringstream ss;
        ss << "Registry<" << typeid(Element).name() << "> with " << _container.size() << " elements";
        return ss.str();
    };

    // получить элемент из регистра
    PYTHON_API std::shared_ptr<Element> getElement(regnum number) {
        if (this->_container.contains(number))
            return this->_container.at(number);
        else
            throw std::invalid_argument(std::format(
                "Такой элемент (<{}> {}) не был зарегистрирован в регистре", typeid(this).name(), number
            ));
    }

    // является ли регистр пустым
    PYTHON_API bool isEmpty() { return _container.empty(); };

    // итератор по элементам регистра
    using iterator = RegistryContainer<Element>::iterator;
    // начало итератора по элементам регистра
    iterator begin() { return iterator(_container.begin()); }
    // конец итератора по элементам регистра
    iterator end() { return iterator(_container.end()); }
    // начало итератора по элементам регистра (константный)
    iterator cbegin() const { return iterator(_container.cbegin()); }
    // конец итератора по элементам регистра (константный)
    iterator cend() const { return iterator(_container.cend()); }
};

// Регистры для материалов и их свойств
//   Предоставляет доступ к элементам, хранящимся в регистре.
using MaterialRegistry = PYTHON_API Registry<Material>;

template <typename Key, typename Values> class Cache {
  private:
    // мьютекс для защиты разделяемых данных (_cache)
    mutable std::shared_mutex _mutex;

  protected:
    // контейнер
    std::unordered_map<Key, Values> _cache;

  public:
    // инициализация всегда пустая
    Cache() : _cache() {};
    ~Cache() = default;

    // Есть ли элемент в кэше
    inline bool has(const Key &key) const {
        std::shared_lock lock(this->_mutex);
        return this->_cache.contains(key);
    }

    // Обновить или добавить элемент в кэш
    inline void update(const Key &key, const Values &value) {
        std::unique_lock lock(this->_mutex);
        this->_cache[key] = value;
    }

    // Получить значение по ключу
    inline const Values &get(const Key &key) const {
        std::shared_lock lock(this->_mutex);
        if (!this->has(key)) {
            throw std::out_of_range("Key {key} not found in cache");
        }
        return this->_cache.at(key);
    }

    // Удалить элемент по ключу
    inline void erase(const Key &key) {
        std::unique_lock lock(this->_mutex);
        this->_cache.erase(key);
    }

    // Удалить все элементы
    inline void clear() {
        std::unique_lock lock(this->_mutex);
        this->_cache.clear();
    }
};

}; // namespace PYTHON_API spindynapy

// ^
// | template classes (registries, caches)
// |
// ================= NEW_BLOCK ===================
// |
// | MACROSES and BINDINGS for PYTHON
// v

#define REGISTRY_TEMPLATE_BINDINGS(cls)                                                                      \
    .def("__str__", &cls::__str__, py::doc("строковое представление для принтинга"))                         \
        .def(                                                                                                \
            "get_element",                                                                                   \
            &cls::getElement,                                                                                \
            py::arg("number"),                                                                               \
            py::return_value_policy::reference,                                                              \
            py::doc("получить элемент из регистра")                                                          \
        )                                                                                                    \
        .def(                                                                                                \
            "is_empty",                                                                                      \
            &cls::isEmpty,                                                                                   \
            "true, если в регистре нет ни одного элемента",                                                  \
            py::doc("является ли регистр пустым")                                                            \
        )

// Функция для привязки регистров к Python.
inline void pyBindRegistries(py::module_ &module) {
    using namespace spindynapy;

    // -------- | REGISTRIES | --------
    py::module_ registries_module = module.def_submodule("registries");

    registries_module.doc() =
        "Регистры - контейнеры, хранящие в себе экземпляры переиспользуемых классов,\n"
        "  которые задают базовые настройки симуляции (например, свойства конкретного материала).\n"
        "\n"
        "У регистров может быть расширенный функционал.\n"
        "\"Легковесы\". Живут всю программу.\n"
        "\n"
        "Кэши - потокобезопасные контейнеры, которые по реализации похожи на регистры,\n"
        "   но предназначены для хранения промежуточных данных, кэш.";

    py::class_<MaterialRegistry, std::shared_ptr<MaterialRegistry>>(registries_module, "MaterialRegistry")
        .def(py::init<>())
        .def(py::init<RegistryContainer<Material>>(), py::arg("container"))
            REGISTRY_TEMPLATE_BINDINGS(MaterialRegistry)
        .doc() = "Регистры для материалов и их свойств\n"
                 "  Предоставляет доступ к элементам, хранящимся в регистре.\n"
                 " БЕЗ ЗАПИСИ => ПОТОКОБЕЗОПАСНЫЙ";
}

#endif // ! __REGISTRIES_HPP__
