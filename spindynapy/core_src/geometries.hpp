#ifndef __GEOMETRIES_BASE_HPP__
#define __GEOMETRIES_BASE_HPP__

/**
 * Базовые интерфейсы геометрии и управлятора макроячеек. ПОТОКОБЕЗОПАСНЫЕ.
 *  И их реализации в (декартовая, ) системах координат.
 *
 *  Геометрия хранит в себе состояние расположения всех моментов и их свойств,
 *    а также предоставляет интерфейс для работы с ними и с макроячейками.
 *  Задаёт внутреннее устройство образца
 *  Геометрия также может поставлять ограниченное количество обсчитываемых параметров.
 */

#include "registries.hpp"
#include "types.hpp"

#include <memory>
#include <object.h>
#include <pybind11/detail/common.h>
#include <pybind11/pybind11.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;

namespace PYTHON_API spindynapy {

// кэш моментов (индекс момента, радиус обрезки) (см. registires.hpp) (ПОТОКОБЕЗОПАСНЫЙ)
using MomentIndexCache = Cache<std::pair<size_t, double>, std::vector<size_t>>;
// кэш макроячеек (индекс макроячейки, радиус обрезки) (см. registries.hpp) со среднеми моментами
// (ПОТОКОБЕЗОПАСНЫЙ)
using MacrocellIndexCache = Cache<std::pair<size_t, double>, std::vector<size_t>>;

// Контейнер для моментов в выбранной системе координат
template <CoordSystemConcept CoordSystem>
using MomentsContainer = PYTHON_API std::vector<std::shared_ptr<typename CoordSystem::Moment>>;

/**
 * Структура макроячейки - это контейнер, который хранит в себе
 *  индексы моментов, которые находятся внутри физической макроячейки - разбиения геометрии
 *
 * Содержит в себе средний момент, который является усреднением всех моментов,
 *  и список индексов моментов, которые находятся внутри макроячейки.
 */
template <CoordSystemConcept CoordSystem> struct PYTHON_API MacroCell {
  public:
    // индексы моментов, находящихся внутри макроячейки
    PYTHON_API std::vector<size_t> moment_indices;

    // средний момент всей макроячейки
    PYTHON_API std::shared_ptr<typename CoordSystem::Moment> avg_moment;

    // конструктор по умолчанию (для создания пустой макроячейки и затем её заполнения)
    MacroCell() : moment_indices(0), avg_moment(nullptr) {};

    // полноценный конструктор
    MacroCell(const std::vector<size_t> &indices, std::shared_ptr<typename CoordSystem::Moment> moment)
        : moment_indices(indices), avg_moment(moment) {};
};

/**
 * Базовый интерфейс управлятора макроячеек.
 *   Предоставляет интерфейс для работы с макроячейками в выбранной системе координат,
 *   включая методы для создания, обновления и получения макроячеек.
 */
template <CoordSystemConcept CoordSystem> class PYTHON_API IMacrocellManager {
  protected:
    // конструктор только для наследников
    IMacrocellManager() = default;

    // очистить набор макроячеек (НЕ ЗАЩИЩЁН МЬЮТЕКСОМ)
    virtual void clearMacrocellsImpl() = 0;

    // создать макроячейки (обнуляя предыдущие) (НЕ ЗАЩИЩЁН МЬЮТЕКСОМ)
    virtual void createMacrocellsImpl() = 0;

    // взять из кэша (НЕ ЗАЩИЩЁН МЬЮТЕКСОМ)
    virtual MomentsContainer<CartesianCoordSystem> getFromCacheImpl(std::pair<size_t, double> cache_key) = 0;

  public:
    // деструктор
    virtual ~IMacrocellManager() = default;

    // создать макроячейки, если они не созданы (или насильно пересоздать)
    PYTHON_API virtual void createMacrocellsIfNotCreated(bool recreate = false) = 0;

    // обновить данные по макроячейкам (средний момент)
    PYTHON_API virtual void updateMacrocells() = 0;

    // получить конейтнер макроячеек
    PYTHON_API virtual const std::vector<MacroCell<CoordSystem>> &getMacrocells() = 0;

    // получить макроячейку по индексу относящегося к ней момента
    PYTHON_API virtual size_t getMacrocellIndexBySpin(size_t spin_index) = 0;

    // обнулить набор макроячеек
    PYTHON_API virtual void clearMacrocells() = 0;

    // получить контейнер моментов из макроячейки по индексу макроячейки (TODO возвращать ссылку из кэша
    // всегда)
    virtual MomentsContainer<CoordSystem> getMomentsFromMacrocells(size_t index, double cutoff_radius) = 0;
};

/**
 * Базовый интерфейс геометрии.
 *   Предоставляет интерфейс для работы с геометрией в выбранной системе координат,
 *   включая методы для получения моментов, их соседей, работы с макроячейками (расширение).
 *
 *   Геометрия хранит в себе состояние расположения всех моментов и их свойств.
 */

template <CoordSystemConcept CoordSystem>
class PYTHON_API IGeometry : virtual public IMacrocellManager<CoordSystem> {
  protected:
    // конструктор только для наследников
    IGeometry() = default;

  public:
    // деструктор
    virtual ~IGeometry() = default;

    // размер геометрии в количестве моментов
    PYTHON_API virtual size_t size() const = 0;

    // возвращает копию моментов из геометрии
    PYTHON_API virtual MomentsContainer<CoordSystem> getMoments() const = 0;

    // предподготовка геометрии (TODO убрать)
    virtual void prepareData() = 0;

    // взять момент по индексу (TODO добавить .at(index))
    PYTHON_API virtual CoordSystem::Moment &operator[](size_t index) = 0;

    // получить контейнер моментов по массиву индексов (TODO добавить .at(vector<index>))
    virtual MomentsContainer<CoordSystem> getFromIndexes(const std::vector<size_t> &indexes) const = 0;

    // получить индексы соседей для момента по индексу в радиусе обрезки (TODO возвращать ссылку из кэша)
    PYTHON_API virtual std::vector<size_t> getNeighbors(size_t index, double cutoff_radius) = 0;

    // клонировать геометрию, создав новый объект (TODO заменить на конструктор копирования)
    virtual std::unique_ptr<IGeometry<CoordSystem>> clone() const = 0;

    // ---- ИТЕРАТОР ----
    // итератор по моментам для обхода
    using iterator = MomentsContainer<CoordSystem>::iterator;
    // начало итератора по моментам
    virtual iterator begin() = 0;
    // конец итератора по моментам
    virtual iterator end() = 0;

    // ---- ПРЕДСТАВЛЕНИЯ ----
    // представить геометрию в виде numpy массива (для сохранения промежуточных результатов)
    PYTHON_API virtual Eigen::MatrixXd asNumpy() const = 0;
    // строковое представление для принтинга (TODO заменить на ostream/stringstream)
    PYTHON_API virtual std::string __str__() const = 0;
};

// ^
// | base template interfaces
// |
// ================= NEW_BLOCK ===================
// |
// | cartesian realization of interfaces (namespace cartesian)
// v

namespace PYTHON_API cartesian {

/**
 * Базовый интерфейс геометрии.
 *   Предоставляет интерфейс для работы с геометрией в выбранной (декартовой) системе координат,
 *   включая методы для получения моментов, их соседей, работы с макроячейками (расширение).
 *
 *   Геометрия хранит в себе состояние расположения всех моментов и их свойств.
 */
using AbstractGeometry = PYTHON_API IGeometry<NamespaceCoordSystem>;

/**
 * Базовый интерфейс управлятора макроячеек.
 *   Предоставляет интерфейс для работы с макроячейками в выбранной (декартовой) системе координат,
 *   включая методы для создания, обновления и получения макроячеек.
 */
using AbstractMacrocellManager = PYTHON_API IMacrocellManager<NamespaceCoordSystem>;

/**
 * Структура макроячейки - это контейнер, который хранит в себе
 *  индексы моментов, которые находятся внутри физической макроячейки - разбиения геометрии
 *
 * Содержит в себе средний момент, который является усреднением всех моментов,
 *  и список индексов моментов, которые находятся внутри макроячейки.
 */
using NamespaceMacroCell = PYTHON_API MacroCell<NamespaceCoordSystem>;

/**
 * Управлятор макроячеек в декартовой системе координат.
 *   Предоставляет интерфейс для работы с макроячейками в выбранной (декартовой) системе координат,
 *   включая методы для создания, обновления и получения макроячеек.
 */
class PYTHON_API MacrocellManager : virtual public AbstractMacrocellManager {
  protected:
    // указатель на контейнер моментов из геометрии
    MomentsContainer<NamespaceCoordSystem> *_moments;
    // кэш соседей (включая свою) - индексов макроячеек с радиусом обрезки
    // (индекс макроячейки, радиус обрезки)
    MacrocellIndexCache _macrocell_index_cache;
    // массив макроячеек (индекс макроячейки -> макроячейка)
    // _macrocells[cell_index] = macrocell meta info
    std::vector<NamespaceMacroCell> _macrocells;
    // карта (индекс спина -> индекс макроячейки (через индекс массива))
    // _spin2cell[spin_index] = cell index for _macrocells
    std::vector<size_t> _spin2cell;

    // мьютекс для защиты разделяемых данных (_spin2cell, _macrocells)
    mutable std::shared_mutex _mutex;

  public:
    // Размер макроячейки (во всех направлениях) в метрах
    PYTHON_API double macrocell_size = 1e-9;

  protected:
    // очистить макроячейки (НЕ ЗАЩИЩЁН МЬЮТЕКСОМ)
    virtual void clearMacrocellsImpl() override {
        this->_macrocells.clear();
        this->_spin2cell.clear();
    }

    // создать макроячейки (обнуляя предыдущие) (НЕ ЗАЩИЩЁН МЬЮТЕКСОМ)
    virtual void createMacrocellsImpl() override {
        // если указатель на моменты не нулевой и массив моментов не пуст
        if (!this->_moments || this->_moments->empty()) {
            this->clearMacrocellsImpl();
            return;
        }

        // размер массива моментов
        size_t num_spins = this->_moments->size();

        // Определение границ геометрии ("что" будем разбивать на ячейки)
        auto min_coords = (*this->_moments)[0]->getCoordinates().asVector();
        auto max_coords = min_coords;

        // найдём покомпонентно минимум и покомпонентно максимум
        for (const auto &moment : *this->_moments) {
            Eigen::Vector3d coords = moment->getCoordinates().asVector();
            min_coords = min_coords.cwiseMin(coords);
            max_coords = max_coords.cwiseMax(coords);
        }

        // добавим дополнительный отступ в половину размера макроячейки
        min_coords -= Eigen::Vector3d::Constant(this->macrocell_size / 2);
        max_coords += Eigen::Vector3d::Constant(this->macrocell_size / 2);

        // размер системы в декартовых измерениях
        Eigen::Vector3d system_size = max_coords - min_coords;

        // Определение размеров сетки (по количеству макроячеек)
        size_t nx = std::max(1UL, static_cast<size_t>(std::ceil(system_size.x() / this->macrocell_size)));
        size_t ny = std::max(1UL, static_cast<size_t>(std::ceil(system_size.y() / this->macrocell_size)));
        size_t nz = std::max(1UL, static_cast<size_t>(std::ceil(system_size.z() / this->macrocell_size)));

        // общее количество макроячеек
        size_t num_cells = nx * ny * nz;

        size_t ix = 0, iy = 0, iz = 0;
        size_t cell_idx = 0;

        this->_spin2cell.reserve(num_spins);

        // Карта для группировки индексов (индекс_ячейки, индексы_спинов)
        std::map<size_t, std::vector<size_t>> spins_by_cell_indicies;

        // Распределение спинов по ячейкам (линеаризация)
        for (size_t spin_idx = 0; spin_idx < num_spins; ++spin_idx) {
            auto &pos = (*this->_moments)[spin_idx]->getCoordinates().asVector();
            ix = std::min(
                nx - 1, static_cast<size_t>(std::floor((pos.x() - min_coords.x()) / this->macrocell_size))
            );
            iy = std::min(
                ny - 1, static_cast<size_t>(std::floor((pos.y() - min_coords.y()) / this->macrocell_size))
            );
            iz = std::min(
                nz - 1, static_cast<size_t>(std::floor((pos.z() - min_coords.z()) / this->macrocell_size))
            );
            // алгоритм линеаризации трёхмерных индексов в ячеечный индекс
            cell_idx = ix + iy * nx + iz * nx * ny;
            // сохраним, что за текущим спином стоит посчитанная ячейка
            this->_spin2cell[spin_idx] = cell_idx;
            // сохраним в карте, что текущий спин относится к cell_idx
            spins_by_cell_indicies[cell_idx].push_back(spin_idx);
        }

        // очищаем макроячейки для предвариельной
        this->_macrocells = {};
        this->_macrocells.reserve(num_cells);

        // для каждой макроячейки k
        for (size_t k = 0; k < num_cells; ++k) {
            // проверяем, что есть хотя бы один спин
            if (spins_by_cell_indicies.contains(k) && !spins_by_cell_indicies[k].empty()) {
                // и записываем (передаём полностью) предсозданную ячейку в финальный список
                MacroCell<NamespaceCoordSystem> new_cell;
                new_cell.moment_indices = spins_by_cell_indicies[k];
                this->_macrocells.push_back(std::move(new_cell));

                // в конечном итоге, записываем за каждым из спинов из карты его ячейку
                auto cell_index = this->_macrocells.size() - 1;
                for (size_t spin_idx : spins_by_cell_indicies[k]) {
                    this->_spin2cell[spin_idx] = cell_index;
                }
            }
        };
    };

    // получить список моментов из индексов макроячеек в радиусе обрезки (НЕ ЗАЩИЩЁН МЬЮТЕКСОМ)
    MomentsContainer<CartesianCoordSystem> MomentsFromMacrocellsImpl(const std::vector<size_t> &macro_indexes
    ) {
        MomentsContainer<CartesianCoordSystem> moments_from_macrocells(macro_indexes.size(), nullptr);
        for (size_t macro_index = 0; macro_index < macro_indexes.size(); ++macro_index) {
            moments_from_macrocells[macro_index] = this->_macrocells[macro_indexes[macro_index]].avg_moment;
        }
        return moments_from_macrocells;
    }

    // проверить кэш (НЕ ЗАЩИЩЁН МЬЮТЕКСОМ)
    virtual MomentsContainer<CartesianCoordSystem> getFromCacheImpl(std::pair<size_t, double> cache_key
    ) override {
        // получаем список "соседних" макроячеек в радиусе обрезки
        auto macro_indexes = this->_macrocell_index_cache.get(cache_key);
        // и возвращаем список моментов из индексов макроячеек в радиусе обрезки
        return this->MomentsFromMacrocellsImpl(macro_indexes);
    };

  public:
    // деструктор
    virtual ~MacrocellManager() = default;

    // полноценный конструктор (принимает указатель на контейнер моментов и размер макроячейки для разбиения)
    MacrocellManager(MomentsContainer<NamespaceCoordSystem> *moments, double macrocell_size)
        : _moments(moments), macrocell_size(macrocell_size) {
        if (macrocell_size <= 0) {
            throw std::invalid_argument("Macrocell size must be positive");
        }
    };

    // получить конейтнер макроячейк
    PYTHON_API virtual const std::vector<NamespaceMacroCell> &getMacrocells() override {
        std::shared_lock lock(_mutex);
        if (this->_macrocells.empty()) {
            throw std::runtime_error("Macrocells are not created. Call prepareData() first.");
        }
        return this->_macrocells;
    };

    // получить макроячейку по индексу относящегося к ней момента (рейзит ошибку, если не совпадают размеры)
    PYTHON_API virtual size_t getMacrocellIndexBySpin(size_t spin_index) override {
        std::shared_lock lock(_mutex);
        if (spin_index >= this->_spin2cell.size()) {
            throw std::out_of_range("Spin index out of range (для него нет макроячейки...)");
        }
        return this->_spin2cell[spin_index];
    };

    // создать макроячейки, если они не созданы (или насильно пересоздать)
    PYTHON_API virtual void createMacrocellsIfNotCreated(bool recreate = false) override {
        std::unique_lock lock(_mutex);
        if (this->_macrocells.empty() || recreate) {
            this->createMacrocellsImpl(); // создать макроячейки, поделив спины
        }
    };

    // получить контейнер моментов из макроячейки по индексу макроячейки
    virtual MomentsContainer<CartesianCoordSystem> getMomentsFromMacrocells(
        size_t index, double cutoff_radius
    ) override {
        std::shared_lock read_cache_lock(_mutex);
        // проверка на нуль
        if (cutoff_radius <= 1e-15)
            return {}; // типичные размеры: 1e-10

        // ---- проверяем кэш ----
        auto cache_key = std::make_pair(index, cutoff_radius);
        if (this->_macrocell_index_cache.has(cache_key)) {
            return this->getFromCacheImpl(cache_key);
        }

        read_cache_lock.unlock(); // разблокируем мьютекс, чтобы не блокировать другие потоки

        // и заблокируем его на запись (другим запрещаем лезть сюда)
        std::unique_lock write_lock(_mutex);

        // ---- проверяем кэш (могли записать за время) ----
        if (this->_macrocell_index_cache.has(cache_key)) {
            return this->getFromCacheImpl(cache_key);
        }

        // ---- иначе честно считаем ----

        if (!this->_moments || this->_moments->empty()) {
            return {};
        }

        std::vector<size_t> neighbors;

        // относительного кого считаем
        const auto &target_coords = (*this->_moments)[index]->getCoordinates();

        // если расстояние от момента до центра макроячейки меньше радиуса обрезки,
        // то сохраняем в "соседей" индекс ячейки
        for (size_t i = 0; i < _macrocells.size(); ++i) {
            const double distance =
                target_coords.getDistanceFrom(_macrocells[i].avg_moment->getCoordinates());
            if (distance <= cutoff_radius) {
                neighbors.push_back(i);
            }
        }

        // ---- сохраняем в кэш ----
        this->_macrocell_index_cache.update(cache_key, neighbors);

        // и возвращаем список моментов из индексов макроячеек в радиусе обрезки
        MomentsContainer<CartesianCoordSystem> moments_from_macrocells(neighbors.size(), nullptr);
        for (size_t macro_index = 0; macro_index < neighbors.size(); ++macro_index) {
            moments_from_macrocells[macro_index] = this->_macrocells[neighbors[macro_index]].avg_moment;
        }
        return moments_from_macrocells;
    };

    // обновить данные по макроячейкам (средний момент)
    PYTHON_API virtual void updateMacrocells() override {
        // блокировка мьютекса
        std::unique_lock lock(_mutex);

        if (!this->_moments || this->_moments->empty()) {
            return;
        }

        // Итерируем по существующим ячейкам
        for (auto &cell : this->_macrocells) {
            size_t spin_count = cell.moment_indices.size();
            if (spin_count == 0)
                continue; // Пропускаем пустые ячейки (их быть так-то не должно)

            // предсоздадим вектора и указатели для среднего момента
            Eigen::Vector3d avg_coordinates_vector(0, 0, 0);
            Eigen::Vector3d avg_direction_vector(0, 0, 0);

            // для определения превосходящего материала (по количеству спинов)
            std::map<std::shared_ptr<Material>, size_t> material_counts;
            std::shared_ptr<Moment> moment = nullptr;
            std::shared_ptr<Material> material = nullptr;

            // "общая" намагниченность насыщения (для определения веса)
            double total_weight = 0.0;

            // обход по всем моментам в ячейке
            for (size_t spin_idx : cell.moment_indices) {
                moment = (*this->_moments)[spin_idx];
                avg_coordinates_vector += moment->getCoordinates().asVector();
                avg_direction_vector += moment->getDirection().asVector();
                if ((material = moment->getMaterialSharedPtr())) {
                    material_counts[material]++;
                    total_weight += material->atomic_magnetic_saturation_magnetization;
                }
            }

            // усреднение на количество спинов
            avg_coordinates_vector /= static_cast<double>(spin_count);
            avg_direction_vector /= static_cast<double>(spin_count);

            // Поиск наиболее часто встречающегося материала по ранее созданной карте
            std::shared_ptr<Material> predominant_material = nullptr;
            size_t max_count = 0;
            for (const auto &[encountered_material, em_count] : material_counts) {
                if (em_count > max_count) {
                    max_count = em_count;
                    predominant_material = encountered_material;
                }
            }

            // мера организованности спинов (если все параллельные)
            auto magnetic_organization_measue =
                avg_direction_vector.norm(); // в интервале [0, 1], где 1 - полностью параллельные

            // наконец, обновляем момент в ячейке
            cell.avg_moment =
                std::make_shared<Moment>(avg_coordinates_vector, avg_direction_vector, predominant_material);
            // также определим "вес" для момента макроячейки (в простейшем случае = количество спинов в
            // ячейке)
            cell.avg_moment->amplitude = magnetic_organization_measue * total_weight /
                                         predominant_material->atomic_magnetic_saturation_magnetization;
        }
    };

    // обнулить набор макроячеек
    PYTHON_API virtual void clearMacrocells() override {
        std::unique_lock lock(_mutex);
        this->clearMacrocellsImpl();
    };
};

/**
 * Геометрия, реализуемая в декартовой системе координат.
 *   Предоставляет интерфейс для работы с геометрией в выбранной (декартовой) системе координат,
 *   включая методы для получения моментов, их соседей, работы с макроячейками (расширение).
 *
 *   Геометрия хранит в себе состояние расположения всех моментов и их свойств.
 */
class Geometry : public AbstractGeometry, public MacrocellManager {
  protected:
    // контейнер моментов в декартовой системе координат
    MomentsContainer<NamespaceCoordSystem> _moments;
    // кэш для соседей-индексов в радиусе обрезки (индекс момента, радиус обрезки) = индексы соседей
    MomentIndexCache _moment_index_cache;

    // мьютекс для защиты разделяемых данных
    mutable std::shared_mutex _mutex;

  public:
    // конструктор изнутри системы
    PYTHON_API Geometry(const MomentsContainer<NamespaceCoordSystem> &moments, double macrocell_size = 1e-9)
        : MacrocellManager(nullptr, macrocell_size), _moments(moments) {
        MacrocellManager::_moments = &this->_moments; // заполнить указатель на массив моментов
    };
    // конструктор изнутри системы для rvalue
    PYTHON_API Geometry(MomentsContainer<NamespaceCoordSystem> &&moments, double macrocell_size = 1e-9)
        : MacrocellManager(nullptr, macrocell_size), _moments(std::move(moments)) {
        MacrocellManager::_moments = &this->_moments; // заполнить указатель на массив моментов
    };
    // главный конструктор: из numpy массива; принимает регистр для проверки материалов, размер макроячейки
    PYTHON_API Geometry(
        const Eigen::MatrixXd &moments, MaterialRegistry &material_registry, double macrocell_size = 1e-9
    )
        : MacrocellManager(nullptr, macrocell_size) {
        //
        if (moments.cols() < 7) {
            throw std::invalid_argument("Expected 7 columns: [x, y, z, sx, sy, sz, material]");
        };
        this->_moments.reserve(static_cast<size_t>(moments.rows()));
        for (Eigen::Index i = 0; i < moments.rows(); ++i) {
            this->_moments.emplace_back(std::make_shared<Moment>(
                Coordinates(moments(i, 0), moments(i, 1), moments(i, 2)),
                Direction(moments(i, 3), moments(i, 4), moments(i, 5)),
                material_registry.getElement(static_cast<regnum>(moments(i, 6)))
            ));
        };
        MacrocellManager::_moments = &this->_moments; // заполнить указатель на массив моментов
    };

    // размер геометрии в количестве моментов
    PYTHON_API virtual size_t size() const override {
        std::shared_lock lock(_mutex);
        return this->_moments.size();
    }

    // возвращает копию моментов из геометрии
    PYTHON_API virtual MomentsContainer<NamespaceCoordSystem> getMoments() const override {
        std::shared_lock lock(_mutex);
        return this->_moments;
    };

    // предподготовка геометрии
    virtual void prepareData() override {
        std::unique_lock lock(_mutex); // мутирует _moments по указателю из менеджера
        this->createMacrocellsIfNotCreated(true); // создать макроячейки, поделив спины
        this->updateMacrocells();                 // посчитать их в первый раз
    };

    // взять момент по индексу
    PYTHON_API virtual Moment &operator[](size_t index) override {
        std::shared_lock lock(_mutex); // TODO может быть медленно!!!
        return *this->_moments[index];
    };

    // получить контейнер моментов по массиву индексов
    virtual MomentsContainer<NamespaceCoordSystem> getFromIndexes(const std::vector<size_t> &indexes
    ) const override {
        std::shared_lock lock(_mutex);
        MomentsContainer<NamespaceCoordSystem> result(indexes.size(), nullptr);
        for (size_t i = 0; i < indexes.size(); ++i) {
            result[i] = this->_moments[indexes[i]];
        }
        return result;
    };

    // получить индексы соседей для момента по индексу в радиусе обрезки
    PYTHON_API virtual std::vector<size_t> getNeighbors(size_t index, double cutoff_radius) override {
        std::shared_lock lock(_mutex);
        // проверка на нуль
        if (cutoff_radius <= 1e-15)
            return {}; // типичные размеры: 1e-10

        // ---- проверяем кэш ----
        auto cache_key = std::make_pair(index, cutoff_radius);
        if (this->_moment_index_cache.has(cache_key)) {
            return this->_moment_index_cache.get(cache_key);
        }

        // ---- иначе честно считаем ----
        std::vector<size_t> neighbors;
        const auto &target_coords = this->_moments[index]->getCoordinates();

        for (size_t i = 0; i < this->_moments.size(); ++i) {
            if (i == index)
                continue; // себя не считаем
            const double distance_sq = target_coords.getDistanceFrom(_moments[i]->getCoordinates());
            if (distance_sq <= cutoff_radius) {
                neighbors.push_back(i);
            }
        }

        // ---- сохраняем в кэш ----
        this->_moment_index_cache.update(cache_key, neighbors);

        // возвращаем индексы соседей
        return neighbors;
    };

    // ---- ИТЕРАТОР ----
    // итератор по моментам для обхода
    using iterator = MomentsContainer<NamespaceCoordSystem>::iterator;
    // начало итератора по моментам
    virtual iterator begin() override {
        std::shared_lock lock(_mutex);
        return this->_moments.begin();
    }
    // конец итератора по моментам
    virtual iterator end() override {
        std::shared_lock lock(_mutex);
        return this->_moments.end();
    }

    // клонировать геометрию, создав новый объект
    virtual std::unique_ptr<IGeometry<NamespaceCoordSystem>> clone() const override {
        std::shared_lock lock(_mutex);
        MomentsContainer<NamespaceCoordSystem> cloned_moments;
        cloned_moments.reserve(this->_moments.size());
        for (const auto &moment : this->_moments) {
            cloned_moments.push_back(moment->clone());
        }
        auto cloned = std::make_unique<Geometry>(cloned_moments, this->macrocell_size);
        cloned->prepareData();
        return cloned;
    }

    // ---- ПРЕДСТАВЛЕНИЯ ----

    // представить геометрию в виде numpy массива (для сохранения промежуточных результатов)
    PYTHON_API virtual Eigen::MatrixXd asNumpy() const override {
        std::shared_lock lock(_mutex);
        auto moments_len = this->_moments.size();
        Eigen::MatrixXd numpy_matrix(moments_len, 7);
        for (size_t i = 0; i < moments_len; ++i) {
            numpy_matrix(i, 0) =
                this->_moments[i]->getCoordinates().asVector().x(); // NOLINT(*-sign-conversion)
            numpy_matrix(i, 1) =
                this->_moments[i]->getCoordinates().asVector().y(); // NOLINT(*-sign-conversion)
            numpy_matrix(i, 2) =
                this->_moments[i]->getCoordinates().asVector().z(); // NOLINT(*-sign-conversion)
            numpy_matrix(i, 3) =
                this->_moments[i]->getDirection().asVector().x(); // NOLINT(*-sign-conversion)
            numpy_matrix(i, 4) =
                this->_moments[i]->getDirection().asVector().y(); // NOLINT(*-sign-conversion)
            numpy_matrix(i, 5) =
                this->_moments[i]->getDirection().asVector().z();              // NOLINT(*-sign-conversion)
            numpy_matrix(i, 6) = this->_moments[i]->getMaterial().getNumber(); // NOLINT(*-sign-conversion)
        }
        return numpy_matrix;
    };

    // строковое представление для принтинга
    PYTHON_API virtual std::string __str__() const override {
        std::shared_lock lock(_mutex);
        std::stringstream ss;
        for (const auto &elem : this->_moments) {
            ss << "\n" << elem->__str__();
        }
        return ss.str();
    };
};

}; // namespace PYTHON_API cartesian

}; // namespace PYTHON_API spindynapy

// ^
// | cartesian realization of interfaces (namespace cartesian)
// |
// ================= NEW_BLOCK ===================
// |
// | MACROSES and GEOMETRIES BINDINGS for PYTHON
// v

#define MACROCELL_TEMPLATE_BINDINGS(cls)                                                                     \
    .def_readonly(                                                                                           \
        "moment_indices", &cls::moment_indices, py::doc("индексы моментов, находящихся внутри макроячейки")  \
    )                                                                                                        \
        .def_property_readonly(                                                                              \
            "avg_moment",                                                                                    \
            [](const cls &self) { return self.avg_moment.get(); },                                           \
            py::return_value_policy::reference_internal,                                                     \
            py::doc("средний момент всей макроячейки")                                                       \
        )

#define ABSTRACTMACROCELLMANAGER_TEMPLATE_BINDINGS(cls)                                                      \
    .def("get_macrocells", &cls::getMacrocells, py::doc("получить конейтнер макроячеек"))                    \
        .def(                                                                                                \
            "get_macrocell_index_by_spin",                                                                   \
            &cls::getMacrocellIndexBySpin,                                                                   \
            py::doc("получить макроячейку по индексу относящегося к ней момента")                            \
        )                                                                                                    \
        .def(                                                                                                \
            "update_macrocells",                                                                             \
            &cls::updateMacrocells,                                                                          \
            py::doc("обновить данные по макроячейкам (средний момент)")                                      \
        )                                                                                                    \
        .def(                                                                                                \
            "create_macrocells_if_not_created",                                                              \
            &cls::createMacrocellsIfNotCreated,                                                              \
            py::arg("recreate") = false,                                                                     \
            py::doc("создать макроячейки, если они не созданы (или насильно пересоздать)")                   \
        )                                                                                                    \
        .def("clear_macrocells", &cls::clearMacrocells, py::doc("обнулить набор макроячеек"))

#define ABSTRACTGEOMETRY_TEMPLATE_BINDINGS(cls)                                                              \
    .def("__str__", &cls::__str__, py::doc("строковое представление для принтинга"))                         \
        .def("__len__", &cls::size, py::doc("размер геометрии в количестве моментов"))                       \
        .def(                                                                                                \
            "get_neighbors",                                                                                 \
            &cls::getNeighbors,                                                                              \
            py::arg("index"),                                                                                \
            py::arg("cutoff_radius"),                                                                        \
            py::doc("получить индексы соседей для момента по индексу в радиусе обрезки")                     \
        )                                                                                                    \
        .def(                                                                                                \
            "as_numpy",                                                                                      \
            &cls::asNumpy,                                                                                   \
            py::doc("представить геометрию в виде numpy массива (для сохранения промежуточных результатов)") \
        )                                                                                                    \
        .def(                                                                                                \
            "__getitem__",                                                                                   \
            &cls::operator[],                                                                                \
            py::arg("index"),                                                                                \
            py::return_value_policy::reference,                                                              \
            py::doc("взять момент по индексу")                                                               \
        )                                                                                                    \
        .def("get_moments", &cls::getMoments, py::doc("получить контейнер моментов из геометрии"))

/**
 * Функция для привязки геометрий к Python.
 *   Создаёт модуль geometries и подмодули систем координат,
 */
inline void pyBindGeometries(py::module_ &module) {
    using namespace spindynapy;

    // -------- | GEOMETRIES | --------
    py::module_ geometries_module = module.def_submodule("geometries");

    geometries_module.doc() =
        "Базовые интерфейсы геометрии и управлятора макроячеек.\n"
        " И их реализации в (декартовая, ) системах координат.\n"
        " Геометрия хранит в себе состояние расположения всех моментов и их свойств,\n"
        "   а также предоставляет интерфейс для работы с ними и с макроячейками.\n"
        " Задаёт внутреннее устройство образца\n"
        " Геометрия также может поставлять ограниченное количество обсчитываемых параметров.";

    // -------- | CARTESIAN GEOMETRIES | --------
    py::module_ cartesian = geometries_module.def_submodule("cartesian");

    cartesian.doc() = "Геометрия в декартовой системе координат.\n"
                      " Задаётся в виде массива моментов, которые хранятся в декартовой системе координат.\n"
                      " Реализует интерфейс геометрии и управлятора макроячеек.\n"
                      " Геометрия хранит в себе состояние расположения всех моментов и их свойств,\n"
                      "   а также предоставляет интерфейс для работы с ними и с макроячейками.\n"
                      " Задаёт внутреннее устройство образца\n"
                      " Геометрия также может поставлять ограниченное количество обсчитываемых параметров.";

    using cartesian::AbstractGeometry;
    using cartesian::AbstractMacrocellManager;
    using cartesian::Geometry;
    using cartesian::MacrocellManager;
    using cartesian::Moment;
    using cartesian::NamespaceMacroCell;

    py::class_<NamespaceMacroCell, std::shared_ptr<NamespaceMacroCell>>(cartesian, "MacroCell")
        MACROCELL_TEMPLATE_BINDINGS(NamespaceMacroCell)
            .doc() =
        "Структура макроячейки - это контейнер, который хранит в себе\n"
        "   индексы моментов, которые находятся внутри физической макроячейки - разбиения геометрии\n"
        "\n"
        "Содержит в себе средний момент, который является усреднением всех моментов,\n"
        "  и список индексов моментов, которые находятся внутри макроячейки.";

    py::class_<AbstractMacrocellManager, std::shared_ptr<AbstractMacrocellManager>>(
        cartesian, "AbstractMacrocellManager"
    ) ABSTRACTMACROCELLMANAGER_TEMPLATE_BINDINGS(AbstractMacrocellManager)
        .doc() =
        "Базовый интерфейс управлятора макроячеек.\n"
        "  Предоставляет интерфейс для работы с макроячейками в выбранной (декартовой) системе координат,\n"
        "  включая методы для создания, обновления и получения макроячеек.\n";

    py::class_<AbstractGeometry, AbstractMacrocellManager, std::shared_ptr<AbstractGeometry>>(
        cartesian, "AbstractGeometry"
    ) ABSTRACTGEOMETRY_TEMPLATE_BINDINGS(AbstractGeometry)
        .doc() =
        "Базовый интерфейс геометрии.\n"
        "  Предоставляет интерфейс для работы с геометрией в выбранной (декартовой) системе координат,\n"
        "  включая методы для получения моментов, их соседей, работы с макроячейками (расширение).\n"
        "\n"
        "  Геометрия хранит в себе состояние расположения всех моментов и их свойств.";

    py::class_<MacrocellManager, AbstractMacrocellManager, std::shared_ptr<MacrocellManager>>(
        cartesian, "MacrocellManager"
    )
        .def_readonly(
            "macrocell_size",
            &MacrocellManager::macrocell_size,
            py::doc("Размер макроячейки (во всех направлениях) в метрах")
        )
        .doc() =
        "Управлятор макроячеек в декартовой системе координат.\n"
        "  Предоставляет интерфейс для работы с макроячейками в выбранной (декартовой) системе координат,\n"
        "  включая методы для создания, обновления и получения макроячеек."
        " (потокобезопасный)";

    py::class_<Geometry, AbstractGeometry, MacrocellManager, std::shared_ptr<Geometry>>(cartesian, "Geometry")
        .def(
            py::init<const MomentsContainer<CartesianCoordSystem> &, double>(),
            py::arg("moments"),
            py::arg("macrocell_size") = 1e-9,
            py::doc("Создаёт геометрию из контейнера моментов (rvalue или lvalue ссылки).\n"
                    "  Принимает размер макроячейки.")
        )
        .def(
            py::init<const Eigen::MatrixXd &, MaterialRegistry &, double>(),
            py::arg("moments"),
            py::arg("material_registry"),
            py::arg("macrocell_size") = 1e-9,
            py::doc("Создаёт геометрию из numpy массива моментов.\n"
                    "  Принимает регистр для проверки материалов и размер макроячейки.")
        )
        .doc() =
        "Геометрия, реализуемая в декартовой системе координат.\n"
        "  Предоставляет интерфейс для работы с геометрией в выбранной (декартовой) системе координат,\n"
        "  включая методы для получения моментов, их соседей, работы с макроячейками (расширение).\n"
        "\n"
        "  Геометрия хранит в себе состояние расположения всех моментов и их свойств."
        " (потокобезопасный)";
}

#endif // ! __GEOMETRIES_BASE_HPP__
