#ifndef __SPINDYNAPY_REGISTRIES_HPP__
#define __SPINDYNAPY_REGISTRIES_HPP__

/**
 * @file   registries.hpp
 * @brief  Read-only registries and thread-safe caches for reusable simulation data.
 *
 * Containers in this header fall into two categories:
 *
 * - Registry – application-lifetime, immutable map keyed by `regnum`
 *      (e.g. predefined materials).  Once populated a registry is never mutated,
 *      so concurrent access is lock-free.
 *
 * - Cache – mutable map protected by `std::shared_mutex`.  Intended for
 *      intermediate data whose loss is acceptable.
 *
 * Exposed entities:
 *   - Template `Registry<Element>` immutable container.
 *   - Alias `MaterialRegistry` specialisation for `Material`.
 *   - Template `Cache<Key, Value>` shared-mutex-protected map.
 *   - `pyBindRegistries()`exports the sub-module `core.registries`.
 *
 * @note  Registries are deliberately lightweight; they must not own heavy
 *           OS resources.  Caches may be cleared at any time without harm.
 *
 * @copyright 2025 SpinDynaPy
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

// ===========================================================================
//  Hash helpers
// ===========================================================================

/**
 * @brief Combine an existing hash `seed` with a new value `key`.
 *
 * Boost-style mix that lowers collision probability for composite keys.
 *
 * @tparam T any type hashable by `std::hash<T>`.
 * @param seed hash accumulator (modified in place).
 * @param key value to mix in.
 * @returns void – mutates `seed`.
 */
template <typename T> void hash_combine(size_t &seed, T const &key) {
    hash<T> hasher;
    seed ^= hasher(key) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/**
 * @brief Order-sensitive hash for `std::pair<T1,T2>`.
 *
 * Guarantees `hash({a,b}) != hash({b,a})` when `a != b`.
 *
 * @tparam T1 first element type.
 * @tparam T2 second element type.
 */
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

// ---------------------------------------------------------------------------
//  Type aliases
// ---------------------------------------------------------------------------

/**
 * @brief Internal storage type for a registry.
 *
 * @details Alias for the map ``{regnum → std::shared_ptr<Element>}``
 *
 * @tparam Element element type stored in the registry.
 */
template <typename Element>
using RegistryContainer = PYTHON_API std::unordered_map<regnum, std::shared_ptr<Element>>;

// ===========================================================================
//  Registry
// ===========================================================================

/**
 * @class Registry
 * @brief Immutable, thread-safe container keyed by `regnum`.
 *
 * @tparam Element stored element type.
 *
 * The internal map is initialised in the constructor and never mutated;
 *   hence look-ups require no locks.
 * Attempting to access a missing key throws `std::invalid_argument`.
 */
template <typename Element> class PYTHON_API Registry {
  protected:
    /** Internal storage for elements. */
    RegistryContainer<Element> _container;

  public:
    /** @brief Construct an empty registry. */
    PYTHON_API Registry() : _container() {};
    /**
     * @brief Construct from a pre-filled container.
     *
     * @param container ready map ``{regnum: Element}``.
     * @returns void – stores the map by move/copy.
     */
    PYTHON_API Registry(RegistryContainer<Element> container) : _container(container) {}
    /** @brief Default destructor (no dynamic resources). */
    ~Registry() = default;

    /**
     * @brief Human-readable summary.
     *
     * @returns `std::string` describing element count.
     */
    PYTHON_API std::string __str__() const {
        std::stringstream ss;
        ss << "Registry<" << typeid(Element).name() << "> with " << _container.size() << " elements";
        return ss.str();
    };

    /**
     * @brief Retrieve element by registry index.
     *
     * @param number registry id (`regnum`).
     * @returns `shared_ptr<Element>` (never null).
     * @throws `std::invalid_argument` if key absent.
     */
    PYTHON_API std::shared_ptr<Element> getElement(regnum number) {
        if (this->_container.contains(number))
            return this->_container.at(number);
        else
            throw std::invalid_argument(std::format(
                "Такой элемент (<{}> {}) не был зарегистрирован в регистре", typeid(this).name(), number
            ));
    }

    /**
     *  @brief Check if the registry is empty.
     *
     *  @returns `true` if registry contains zero elements.
     */
    PYTHON_API bool isEmpty() { return _container.empty(); };

    /// @brief Начало (неконстантное).
    using iterator = RegistryContainer<Element>::iterator;
    /// @brief Начало (неконстантное).
    iterator begin() { return iterator(_container.begin()); }
    /// @brief Конец (неконстантное).
    iterator end() { return iterator(_container.end()); }
    /// @brief Начало (константное).
    iterator cbegin() const { return iterator(_container.cbegin()); }
    /// @brief Конец (константное).
    iterator cend() const { return iterator(_container.cend()); }
};

/**
 * @brief Convenience alias for a registry of predefined `Material` instances.
 *
 * Populated once at start-up; subsequent access is read-only and therefore lock-free.
 */
using MaterialRegistry = PYTHON_API Registry<Material>;

// ===========================================================================
//  Cache
// ===========================================================================

/**
 * @class   Cache
 * @brief   Mutable map protected by `std::shared_mutex`.
 *
 * @tparam  Key     key type (hashable).
 * @tparam  Value   stored value type.
 *
 * *Readers* (`has`, `get`) acquire a shared lock.
 * *Writers* (`update`, `erase`, `clear`) acquire an exclusive lock.
 */
template <typename Key, typename Values> class Cache {
  private:
    /** @brief Guards all access to `_cache`. */
    mutable std::shared_mutex _mutex;

  protected:
    /** @brief Key→value map protected by `_mutex`. */
    std::unordered_map<Key, Values> _cache;

  public:
    /** @brief Construct an empty registry (always). */
    Cache() : _cache() {};
    /** @brief Default destructor (no dynamic resources). */
    ~Cache() = default;

    /**
     * @brief Check presence of *key*.
     *
     * @param key lookup key.
     * @returns `true` if entry exists.
     */
    inline bool has(const Key &key) const {
        std::shared_lock lock(this->_mutex);
        return this->_cache.contains(key);
    }

    /**
     * @brief Insert or overwrite entry.
     *
     * @param key lookup key.
     * @param value data to store.
     * @returns void – mutates cache.
     */
    inline void update(const Key &key, const Values &value) {
        std::unique_lock lock(this->_mutex);
        this->_cache[key] = value;
    }

    /**
     * @brief Retrieve value; throws if absent.
     *
     * @param key lookup key.
     * @returns const reference to stored value.
     * @throws `std::out_of_range` if key absent.
     */
    inline const Values &get(const Key &key) const {
        std::shared_lock lock(this->_mutex);
        if (!this->has(key)) {
            throw std::out_of_range("Key {key} not found in cache");
        }
        return this->_cache.at(key);
    }

    /**
     * @brief Retrieve value or return fallback.
     *
     * @param key lookup key.
     * @param fallback value to return if key absent.
     * @returns stored value or *fallback* (no insertion).
     */
    inline const Values &get(const Key &key, Values _default) const {
        std::shared_lock lock(this->_mutex);
        if (!this->has(key)) {
            return _default;
        }
        return this->_cache.at(key);
    }

    /**
     * @brief Remove entry (silently ignored if absent).
     *
     * @param key lookup key.
     * @returns void – mutates cache.
     */
    inline void erase(const Key &key) {
        std::unique_lock lock(this->_mutex);
        this->_cache.erase(key);
    }

    /**
     * @brief Clear entire cache.
     *
     * @returns void – mutates cache.
     */
    inline void clear() {
        std::unique_lock lock(this->_mutex);
        this->_cache.clear();
    }
};

}; // namespace PYTHON_API spindynapy

/** @brief Inject common bindings into a *Registry* class definition. */
#define REGISTRY_TEMPLATE_BINDINGS(cls)                                                                      \
    .def("__str__", &cls::__str__, py::doc("@brief Text summary (type + size).\n\n@returns ``str``."))       \
        .def(                                                                                                \
            "get_element",                                                                                   \
            &cls::getElement,                                                                                \
            py::arg("number"),                                                                               \
            py::return_value_policy::reference,                                                              \
            py::doc("@brief Retrieve element by registry id.\n"                                              \
                    "\n"                                                                                     \
                    "@param number ``regnum``.\n"                                                            \
                    "@returns stored element.")                                                              \
        )                                                                                                    \
        .def(                                                                                                \
            "is_empty",                                                                                      \
            &cls::isEmpty,                                                                                   \
            "true, если в регистре нет ни одного элемента",                                                  \
            py::doc("@brief Check if the registry is empty.\n\n"                                             \
                    "@returns `true` if registry contains zero elements.")                                   \
        )

/**
 * @brief Bind registries to a Python sub‑module.
 *
 * @param module Parent PyBind11 module (usually the core extension module).
 * @returns void – extends the parent module.
 */
inline void pyBindRegistries(py::module_ &module) {
    using namespace spindynapy;

    // ---------------------------------------------------------------------
    //  Create submodule
    // ---------------------------------------------------------------------

    py::module_ registries_module = module.def_submodule("registries");

    registries_module.doc() =
        "Read-only registries and thread-safe caches for reusable simulation data.\n"
        "\n"
        "Containers in this header fall into two categories:\n"
        "\n"
        "- Registry – application-lifetime, immutable map keyed by `regnum`\n"
        "    (e.g. predefined materials).  Once populated a registry is never mutated,\n"
        "    so concurrent access is lock-free.\n"
        "\n"
        "- Cache – mutable map protected by `std::shared_mutex`.  Intended for\n"
        "    intermediate data whose loss is acceptable.\n"
        "\n"
        "Exposed entities:\n"
        "  - Template `Registry<Element>` immutable container.\n"
        "  - Alias `MaterialRegistry` specialisation for `Material`.\n"
        "  - Template `Cache<Key, Value>` shared-mutex-protected map.\n"
        "  - `pyBindRegistries()`exports the sub-module `core.registries`.\n"
        "\n"
        "@note  Registries are deliberately lightweight; they must not own heavy\n"
        "       OS resources.  Caches may be cleared at any time without harm.\n"
        "\n"
        "@copyright 2025 SpinDynaPy";

    // ---------------------------------------------------------------------
    //  MaterialRegistry binding
    // ---------------------------------------------------------------------

    py::class_<MaterialRegistry, std::shared_ptr<MaterialRegistry>>(registries_module, "MaterialRegistry")
        .def(py::init<>())
        .def(py::init<RegistryContainer<Material>>(), py::arg("container"))
            REGISTRY_TEMPLATE_BINDINGS(MaterialRegistry)
        .doc() = "@class Registry\n"
                 "@brief Immutable, thread-safe container keyed by `regnum`.\n"
                 "       Convenience alias for a registry of predefined `Material` instances.\n"
                 "\n"
                 "@tparam Element stored element type.\n"
                 "\n"
                 "The internal map is initialised in the constructor and never mutated;\n"
                 "  hence look-ups require no locks.\n"
                 "Attempting to access a missing key throws `std::invalid_argument`.";
}

#endif // ! __SPINDYNAPY_REGISTRIES_HPP__
