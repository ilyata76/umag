#ifndef __GEOMETRIES_BASE_HPP__
#define __GEOMETRIES_BASE_HPP__

/**
 * @file   geometries_base.hpp
 * @brief  Geometry interfaces and macrocell management for spin dynamics simulations.
 *
 * This header defines the abstract interfaces and concrete implementations
 *   for spatial geometry and macrocell management used in spin dynamics simulations.
 *
 * Functional overview:
 *
 * - A **Geometry** manages the spatial layout of magnetic moments in a sample.
 *   It provides access to moment positions, directions, and materials.
 *   It also allows neighborhood queries and macrocell-based approximations.
 *
 * - A **MacrocellManager** subdivides the geometry into spatial regions ("macrocells"),
 *   which aggregate moments for use in demagnetization and dipolar field calculations.
 *
 * - Thread safety:
 *   - All geometry implementations (e.g., `Geometry` in Cartesian coordinates)
 *     are thread-safe via shared mutexes.
 *   - Neighbor and macrocell lookups are cached using thread-safe caches.
 *
 * Interfaces:
 * - `IGeometry<CoordSystem>`         – interface for geometry access and manipulation,
 * - `IMacrocellManager<CoordSystem>`– interface for macrocell-based spatial aggregation,
 * - `MacroCell<CoordSystem>`        – container for grouped moment indices and average moment,
 * - `MomentsContainer<CoordSystem>` – alias for list of shared pointers to moments.
 *
 * Concrete (Cartesian):
 * - `Geometry`                       – thread-safe geometry class in Cartesian system,
 * - `MacrocellManager`              – thread-safe manager for Cartesian macrocells,
 * - `MomentIndexCache`              – cache for moment neighbors by radius,
 * - `MacrocellIndexCache`           – cache for macrocell neighbors by radius.
 *
 * Python bindings:
 * - Macros:
 *   - `MACROCELL_TEMPLATE_BINDINGS`
 *   - `MACROCELLMANAGER_TEMPLATE_BINDINGS`
 *   - `GEOMETRY_TEMPLATE_BINDINGS`
 * - Binder function:
 *   - `pyBindGeometries(py::module_ &)`
 *
 * @note All units are SI unless explicitly stated. Coordinates are in meters.
 *       Directions are dimensionless unit vectors. Fields are in [T], energies in [J].
 *
 * @note This header only covers Cartesian coordinates; other systems may extend it in the future.
 *
 * @copyright 2025 SpinDynaPy
 */

#include "logger.hpp"
#include "registries.hpp"
#include "types.hpp"

#include <format>
#include <memory>
#include <object.h>
#include <ostream>
#include <pybind11/detail/common.h>
#include <pybind11/pybind11.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;

namespace PYTHON_API spindynapy {

/**
 * @typedef MomentIndexCache
 * @brief Thread-safe cache for neighbor indices in the geometry.
 *
 * Caches the neighbor indices for a given moment index and cutoff radius:
 *   (moment_index, cutoff_radius) → vector of neighbor indices.
 *
 * @note Ensures efficient reuse of neighborhood queries across time steps or evaluations.
 */
using MomentIndexCache = Cache<std::pair<size_t, double>, std::vector<size_t>>;

/**
 * @typedef MacrocellIndexCache
 * @brief Thread-safe cache for macrocell neighborhood indices.
 *
 * Caches the neighbor macrocells for a given macrocell index and cutoff radius:
 *   (macrocell_index, cutoff_radius) → vector of neighbor macrocell indices.
 *
 * @note Used to accelerate demagnetization and dipole-dipole interactions in macrocell approximation.
 */
using MacrocellIndexCache = Cache<std::pair<size_t, double>, std::vector<size_t>>;

/**
 * @typedef MomentsContainer
 * @brief Container of magnetic moments in a given coordinate system.
 *
 * Represents a collection of shared pointers to `Moment` objects for a specific coordinate system.
 *
 * @tparam CoordSystem The coordinate system defining the geometry (e.g., Cartesian).
 *
 * @note This container is used across geometry, interactions, and simulation stages.
 */
template <CoordSystemConcept CoordSystem>
using MomentsContainer = PYTHON_API std::vector<std::shared_ptr<typename CoordSystem::Moment>>;

// ===========================================================================
//  MacroCell
// ===========================================================================

/**
 * @struct MacroCell
 * @brief Macrocell container structure for grouping magnetic moments.
 *
 * Represents a logical subdivision of the geometry into macrocells. Each macrocell aggregates a subset
 *   of magnetic moments and stores their average moment for use in approximated long-range interactions
 *   such as demagnetization or dipole-dipole coupling.
 *
 * The macrocell contains:
 *   - A list of indices of moments it owns.
 *   - A shared pointer to the average magnetic moment computed from all internal moments.
 *
 * @tparam CoordSystem The coordinate system used for the geometry (e.g., Cartesian).
 *
 * @note The average moment is updated via `updateMacrocells()` and is used in field approximations.
 */
template <CoordSystemConcept CoordSystem> struct PYTHON_API MacroCell {
  public:
    /**
     * @brief Indices of moments contained within this macrocell.
     *
     * Represents the positions (indices) in the global geometry container that belong to this cell.
     */
    PYTHON_API std::vector<size_t> moment_indices;

    /**
     * @brief Average moment of the macrocell.
     *
     * Computed from all internal moments using vector averaging of position and direction.
     * @details Used in demagnetization, dipole-dipole interactions, etc.
     */
    PYTHON_API std::shared_ptr<typename CoordSystem::Moment> avg_moment;

    /**
     * @brief Default constructor for empty macrocell.
     *
     * Initializes an empty macrocell with no members and null average
     *   for use in initialization and subsequent population.
     */
    MacroCell() : moment_indices(0), avg_moment(nullptr) {};

    /**
     * @brief Construct a macrocell from given moment indices and an average moment.
     *
     * @param indices Indices of moments assigned to this macrocell.
     * @param moment  Average moment (typically computed externally).
     */
    MacroCell(const std::vector<size_t> &indices, std::shared_ptr<typename CoordSystem::Moment> moment)
        : moment_indices(indices), avg_moment(moment) {};
};

// ===========================================================================
//  MacroCells Manager
// ===========================================================================

/**
 * @class IMacrocellManager
 * @brief Abstract interface for macrocell management in a geometry.
 *
 * This interface defines operations for managing macrocells in a geometry (moments lists)
 *   for a given coordinate system. Macrocells are logical groups of spins
 *   used to accelerate long-range interactions via spatial aggregation.
 *
 * Responsibilities include:
 *  - Macrocell creation and re-creation.
 *  - Updating average spin states within cells.
 *  - Mapping between spins and their macrocells.
 *  - Retrieving moments within nearby macrocells
 * etc.
 *
 * @tparam CoordSystem The coordinate system used for the geometry (e.g., Cartesian).
 *
 * @note This interface is extended by concrete geometry implementations and provides
 *       internal mechanisms to interact with the macrocell cache.
 *
 * @see MacroCell
 */
template <CoordSystemConcept CoordSystem> class PYTHON_API IMacrocellManager {
  protected:
    /**
     * @brief Protected default constructor.
     * Only accessible to subclasses implementing macrocell logic.
     */
    IMacrocellManager() = default;

    /**
     * @brief Clear all macrocells (unsafe, no locking).
     *
     * Used by internal mechanisms. Should be guarded externally by locks.
     * @warning thread unsafe
     *
     * @returns void - mutates internal macrocell containers.
     */
    virtual void clearMacrocellsImpl() = 0;

    /**
     * @brief Construct macrocells from moment data (unsafe, no locking).
     *
     * Resets and rebuilds macrocell layout. Should be guarded externally by locks.
     * @warning thread unsafe
     *
     * @returns void - mutates internal macrocell containers.
     */
    virtual void createMacrocellsImpl() = 0;

    /**
     * @brief Retrieve cached macrocell-based neighbor moments (unsafe).
     *
     * Returns a container of averaged macrocell moments within the cutoff radius.
     * Assumes the cache key is valid.
     * @warning thread unsafe
     * @param cache_key Pair of {spin index, cutoff radius}.
     *
     * @returns MomentsContainer from macrocell cache.
     */
    virtual MomentsContainer<CartesianCoordSystem> getFromCacheImpl(std::pair<size_t, double> cache_key) = 0;

  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IMacrocellManager() = default;

    /**
     * @brief Create macrocells if not already present, or recreate if requested.
     *
     * @param recreate Force regeneration of macrocell structure (default: false).
     *
     * @returns void - mutates internal state if needed.
     */
    PYTHON_API virtual void createMacrocellsIfNotCreated(bool recreate = false) = 0;

    /**
     * @brief Recalculate average values for each macrocell.
     *
     * Updates properties (average position, direction, material, amplitude, etc.) for each cell.
     *
     * @returns void - mutates macrocell average moment state.
     */
    PYTHON_API virtual void updateMacrocells() = 0;

    /**
     * @brief Access the full macrocell container.
     *
     * @returns Reference to internal vector of macrocells.
     */
    PYTHON_API virtual const std::vector<MacroCell<CoordSystem>> &getMacrocells() = 0;

    /**
     * @brief Get the macrocell index for a given spin index.
     *
     * Maps a spin index to its owning macrocell.
     * @param spin_index Index of spin in global geometry.
     *
     * @returns Index of macrocell containing this spin.
     */
    PYTHON_API virtual size_t getMacrocellIndexBySpin(size_t spin_index) = 0;

    /**
     * @brief Clear macrocell state.
     *
     * @returns void - mutates internal macrocell containers.
     */
    PYTHON_API virtual void clearMacrocells() = 0;

    /**
     * @brief Get linear macrocell size (length of side) in meters.
     *
     * @returns Edge length of macrocell in SI units [m].
     */
    virtual double getMacrocellSize() const = 0;

    /**
     * @brief Get macrocell volume in cubic meters.
     *
     * @returns Volume of a single macrocell in SI units [m³].
     */
    virtual double getMacrocellVolume() const = 0;

    /**
     * @brief Retrieve neighboring macrocell moments within a cutoff radius.
     *
     * Performs spatial filtering of macrocells and returns their average moments.
     * Used in demagnetization/dipolar approximations for example.
     *
     * @param index          Spin index serving as query center.
     * @param cutoff_radius  Search radius [m].
     *
     * @returns Container of shared pointers to neighboring macrocell moments.
     */
    virtual MomentsContainer<CoordSystem> getMomentsFromMacrocells(size_t index, double cutoff_radius) = 0;
};

// ===========================================================================
//  Geometry
// ===========================================================================

/**
 * @class IGeometry
 * @brief Abstract interface for geometry in a given coordinate system.
 *
 * This interface defines operations for accessing and manipulating
 *   the physical layout of magnetic moments in a system. It also extends
 *   macrocell-related functionality from `IMacrocellManager`.
 *
 * Geometry instances own and expose:
 *   - All magnetic moments in the sample.
 *   - Spatial layout and neighbor information.
 *   - Macrocell-based spatial grouping and aggregation.
 *   - Read-only or cloned access to internal state.
 *
 * @tparam CoordSystem The coordinate system in which the geometry is defined (e.g., Cartesian).
 *
 * @note All units are SI.
 */
template <CoordSystemConcept CoordSystem>
class PYTHON_API IGeometry : virtual public IMacrocellManager<CoordSystem> {
  protected:
    /**
     * @brief Protected constructor.
     * Intended for subclass implementation only.
     */
    IGeometry() = default;

  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IGeometry() = default;

    /**
     * @brief Get number of magnetic moments in the geometry.
     *
     * @returns Number of moments (dimension of the system).
     */
    PYTHON_API virtual size_t size() const = 0;

    /**
     * @brief Get a copy of all moments in the geometry.
     *
     * Returns a container of shared pointers to internal moment objects.
     * Can be used for analysis or iteration.
     *
     * @returns Container of all moments in the system.
     */
    PYTHON_API virtual MomentsContainer<CoordSystem> getMoments() const = 0;

    /**
     * @brief Optional geometry preparation step.
     *
     * This may construct internal spatial structures, create macrocells,
     *   or cache neighbor relations. May mutate internal state.
     * @note Will be replaced by proper constructor-based initialization in future.
     *
     * @returns void - mutates geometry and macrocell state.
     */
    virtual void prepare() = 0;

    /**
     * @brief Access a mutable reference to a moment by index.
     *
     * @param index Index of the moment in the geometry container.
     *
     * @returns Reference to mutable moment.
     */
    PYTHON_API virtual CoordSystem::Moment &operator[](size_t index) = 0;

    /**
     * @brief Extract moments from specified indices.
     *
     * @param indexes Vector of moment indices to extract.
     *
     * @returns Container of corresponding moments.
     */
    virtual MomentsContainer<CoordSystem> getFromIndexes(const std::vector<size_t> &indexes) const = 0;

    /**
     * @brief Get neighboring moment indices for a given moment.
     *
     * Performs spatial filtering using Euclidean distance.
     *
     * @param index         Index of the center moment.
     * @param cutoff_radius Cutoff radius for neighbor search [m].
     *
     * @returns Vector of neighbor indices.
     */
    PYTHON_API virtual std::vector<size_t> getNeighbors(size_t index, double cutoff_radius) = 0;

    /**
     * @brief Clone this geometry and return a new instance.
     *
     * Can optionally share internal caches and trigger preparation steps.
     *
     * @param share_cache Whether to reuse internal caches from the original.
     * @param prepare     Whether to invoke `prepare()` after cloning.
     *
     * @returns Unique pointer to new geometry instance.
     */
    virtual std::unique_ptr<IGeometry<CoordSystem>> clone(bool share_cache = false, bool prepare = false)
        const = 0;

    // ---------------------------------------------------------------------
    //  Iterators
    // ---------------------------------------------------------------------

    /**
     * @brief Type alias for iterator over internal moments.
     */
    using iterator = MomentsContainer<CoordSystem>::iterator;

    /**
     * @brief Iterator to the beginning of the moment container.
     *
     * @returns Iterator to first moment.
     */
    virtual iterator begin() = 0;

    /**
     * @brief Iterator to the end of the moment container.
     *
     * @returns Iterator past the last moment.
     */
    virtual iterator end() = 0;

    // ---------------------------------------------------------------------
    //  Representations
    // ---------------------------------------------------------------------

    /**
     * @brief Export geometry as NumPy-like matrix.
     *
     * Each row contains:
     *   coords, direction, material, etc.
     *
     * @returns Matrix of shape (N, xxx).
     */
    PYTHON_API virtual Eigen::MatrixXd asNumpy() const = 0;

    /**
     * @brief Get human-readable string representation.
     *
     * Typically used for debugging or logging.
     *
     * @returns Stringified geometry contents.
     */
    PYTHON_API virtual std::string __str__() const = 0;
};

// ===========================================================================
//  Cartesians
// ===========================================================================

namespace PYTHON_API cartesian {

/**
 * @class AbstractGeometry
 * @brief Abstract interface for geometry in a cartesian coordinate system.
 *
 * This interface defines operations for accessing and manipulating
 *   the physical layout of magnetic moments in a system. It also extends
 *   macrocell-related functionality from `IMacrocellManager`.
 *
 * Geometry instances own and expose:
 *   - All magnetic moments in the sample.
 *   - Spatial layout and neighbor information.
 *   - Macrocell-based spatial grouping and aggregation.
 *   - Read-only or cloned access to internal state.
 *
 * @note All units are SI.
 */
using AbstractGeometry = PYTHON_API IGeometry<NamespaceCoordSystem>;

/**
 * @class AbstractMacrocellManager
 * @brief Abstract interface for macrocell management in a geometry.
 *
 * This interface defines operations for managing macrocells in a geometry (moments lists)
 *   for a given coordinate system. Macrocells are logical groups of spins
 *   used to accelerate long-range interactions via spatial aggregation.
 *
 * Responsibilities include:
 *  - Macrocell creation and re-creation.
 *  - Updating average spin states within cells.
 *  - Mapping between spins and their macrocells.
 *  - Retrieving moments within nearby macrocells
 * etc.
 *
 * @note This interface is extended by concrete geometry implementations and provides
 *       internal mechanisms to interact with the macrocell cache.
 *
 * @see MacroCell
 */
using AbstractMacrocellManager = PYTHON_API IMacrocellManager<NamespaceCoordSystem>;

/**
 * @struct MacroCell
 * @brief Macrocell container structure for grouping magnetic moments.
 *
 * Represents a logical subdivision of the geometry into macrocells. Each macrocell aggregates a subset
 *   of magnetic moments and stores their average moment for use in approximated long-range interactions
 *   such as demagnetization or dipole-dipole coupling.
 *
 * The macrocell contains:
 *   - A list of indices of moments it owns.
 *   - A shared pointer to the average magnetic moment computed from all internal moments.
 *
 * @note The average moment is updated via `updateMacrocells()` and is used in field approximations.
 */
using NamespaceMacroCell = PYTHON_API MacroCell<NamespaceCoordSystem>;

/**
 * @class MacrocellManager
 * @brief Thread-safe manager for spatial macrocell decomposition in Cartesian coordinate systems.
 *
 * Provides creation, caching, and access to macrocells, which group atomic spins into
 *   spatial regions used for approximated long-range interactions such as demagnetization.
 *
 * Internally uses a regular grid aligned with the geometry bounding box, with
 *   fixed macrocell size. Supports cache-aware neighbor lookup and shared cache usage.
 *
 * Responsibilities include:
 *  - Macrocell creation and re-creation.
 *  - Updating average spin states within cells.
 *  - Mapping between spins and their macrocells.
 *  - Retrieving moments within nearby macrocells
 * etc.
 *
 * @note This class assumes a fixed-size cubic macrocell and uses a thread-safe shared cache
 *   for neighbor lookups in the cutoff radius.
 */
class PYTHON_API MacrocellManager : virtual public AbstractMacrocellManager {
  protected:
    /**
     * @brief Pointer to the container of spin moments from the geometry.
     *
     * Used during macrocell creation and update.
     * @note This is a shared pointer to the moments container
     *       to allow shared access across multiple managers and geometry.
     * @note The moments container is expected to be non-empty when creating macrocells.
     */
    MomentsContainer<NamespaceCoordSystem> *_moments;

    /**
     * @brief Cache for macrocell neighbor lookups by cutoff radius.
     *
     * Maps (cell index, cutoff radius) → list of neighbor macrocell indices.
     * @note This is a shared pointer to allow shared access across multiple managers.
     */
    std::shared_ptr<MacrocellIndexCache> _macrocell_index_cache;

    /**
     * @brief Internal storage of macrocells.
     *
     * Each macrocell holds moment indices and an average moment.
     * @note This is a shared pointer to allow shared access across multiple managers.
     */
    std::shared_ptr<std::vector<NamespaceMacroCell>> _macrocells;

    /**
     * @brief Mapping from spin index to owning macrocell index.
     *
     * Used for fast lookup during interaction evaluation.
     * Maps spin indices to their corresponding macrocell indices.
     * @note This is a shared pointer to allow shared access across multiple managers.
     */
    std::shared_ptr<std::vector<size_t>> _spin2cell;

    /**
     * @brief Shared mutex for concurrent access to macrocell state.
     *
     * Protecs internal state during macrocell creation, updates, and access.
     */
    mutable std::shared_mutex _mutex;

  public:
    /**
     * @brief Size of each macrocell in all directions [m].
     *
     * @note Macrocell is a cubic structure, so this defines the edge length.
     * Defines the granularity of spatial decomposition.
     */
    PYTHON_API double macrocell_size = 1e-9;

  protected:
    /**
     * @brief Clear all macrocells and spin-to-cell mapping.
     *
     * Removes all macrocells and internal indexing data.
     * Does not acquire locks – caller must ensure synchronization.
     *
     * @warning Not thread-safe on its own. Intended for internal use.
     *          Should be called under a lock to ensure thread safety.
     *
     * @returns void - mutates internal macrocell containers.
     */
    virtual void clearMacrocellsImpl() override {
        SCOPED_LOG_TIMER_DEBUG("│  ├─ Clearing macrocells");
        this->_macrocells->clear();
        this->_spin2cell->clear();
    }

    /**
     * @brief Build spatial macrocell decomposition from the geometry.
     *
     * Divides spins into regularly spaced cubic regions, builds
     *   average macrocells, and populates indexing structures.
     *
     * @note Overwrites all existing macrocell data.
     *
     * @warning Not thread-safe on its own. Caller must hold exclusive lock.
     *          Should be called under a lock to ensure thread safety.
     *
     * @throws std::runtime_error if geometry pointer is null or empty.
     *
     * @returns void - mutates internal macrocell containers and mappings.
     */
    virtual void createMacrocellsImpl() override {
        SCOPED_LOG_TIMER_DEBUG("│  ├─ Creating macrocells");
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

        // запись по индексу будет
        this->_spin2cell->resize(num_spins);

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
            this->_spin2cell->at(spin_idx) = cell_idx;
            // сохраним в карте, что текущий спин относится к cell_idx
            spins_by_cell_indicies[cell_idx].push_back(spin_idx);
        }

        // очищаем макроячейки для предвариельной
        this->_macrocells = std::make_shared<std::vector<NamespaceMacroCell>>();
        this->_macrocells->reserve(num_cells);

        // для каждой макроячейки k
        for (size_t k = 0; k < num_cells; ++k) {
            // проверяем, что есть хотя бы один спин
            if (spins_by_cell_indicies.contains(k) && !spins_by_cell_indicies[k].empty()) {
                // и записываем (передаём полностью) предсозданную ячейку в финальный список
                MacroCell<NamespaceCoordSystem> new_cell;
                new_cell.moment_indices = spins_by_cell_indicies[k];
                this->_macrocells->push_back(std::move(new_cell));

                // в конечном итоге, записываем за каждым из спинов из карты его ячейку
                auto cell_index = this->_macrocells->size() - 1;
                for (size_t spin_idx : spins_by_cell_indicies[k]) {
                    this->_spin2cell->at(spin_idx) = cell_index;
                }
            }
        };
    };

    /**
     * @brief Retrieve average moments from a set of macrocell indices.
     *
     * For each macrocell index provided, this method extracts the associated
     *   average moment pointer and returns them as a flat container.
     *
     * This utility is used when constructing a set of neighbor macrocells for
     *   interactions such as dipole-dipole or demagnetization.
     *
     * @param macro_indexes A list of macrocell indices from which to fetch average moments.
     *
     * @note This method is not thread-safe and assumes caller holds appropriate lock.
     *       Should called under a lock to ensure thread safety.
     *
     * @returns Vector of shared pointers to average magnetic moments.
     */
    MomentsContainer<CartesianCoordSystem> MomentsFromMacrocellsImpl(const std::vector<size_t> &macro_indexes
    ) {
        MomentsContainer<CartesianCoordSystem> moments_from_macrocells(macro_indexes.size(), nullptr);
        for (size_t macro_index = 0; macro_index < macro_indexes.size(); ++macro_index) {
            moments_from_macrocells[macro_index] =
                this->_macrocells->at(macro_indexes[macro_index]).avg_moment;
        }
        return moments_from_macrocells;
    }

    /**
     * @brief Lookup cached average moments from macrocells near a spin.
     *
     * Queries internal neighbor cache using a key of (spin index, cutoff radius),
     *   and returns the corresponding precomputed average moments.
     *
     * @param cache_key Pair {spin index, cutoff radius}.
     *
     * @note Does not perform any synchronization. Caller must hold lock.
     *
     * @returns Container of average moment pointers for neighbor macrocells.
     */
    virtual MomentsContainer<CartesianCoordSystem> getFromCacheImpl(std::pair<size_t, double> cache_key
    ) override {
        // получаем список "соседних" макроячеек в радиусе обрезки
        auto macro_indexes = this->_macrocell_index_cache->get(cache_key);
        // и возвращаем список моментов из индексов макроячеек в радиусе обрезки
        return this->MomentsFromMacrocellsImpl(macro_indexes);
    };

  public:
    /**
     * @brief Destructor (no dynamic resources to release).
     */
    virtual ~MacrocellManager() = default;

    /**
     * @brief Construct a macrocell manager from an external moment container.
     *
     * @param moments         Pointer to the container of moments.
     * @param macrocell_size  Spatial resolution for cubic macrocell grid [m].
     *
     * @throws std::invalid_argument.
     */
    MacrocellManager(MomentsContainer<NamespaceCoordSystem> *moments, double macrocell_size)
        : _moments(moments),
          _macrocell_index_cache(std::make_shared<MacrocellIndexCache>()),
          _macrocells(std::make_shared<std::vector<NamespaceMacroCell>>()),
          _spin2cell(std::make_shared<std::vector<size_t>>()),
          macrocell_size(macrocell_size) {
        if (macrocell_size <= 0) {
            throw std::invalid_argument("Macrocell size must be positive");
        }
    };

    /**
     * @brief Access the full macrocell container.
     *
     * @returns Reference to internal vector of macrocells.
     */
    PYTHON_API virtual const std::vector<NamespaceMacroCell> &getMacrocells() override {
        std::shared_lock lock(_mutex);
        if (this->_macrocells->empty()) {
            throw std::runtime_error("Macrocells are not created. Call prepareData() first.");
        }
        return *this->_macrocells.get();
    };

    /**
     * @brief Get the macrocell index for a given spin index.
     *
     * Maps a spin index to its owning macrocell.
     * @param spin_index Index of spin in global geometry.
     *
     * @returns Index of macrocell containing this spin.
     */
    PYTHON_API virtual size_t getMacrocellIndexBySpin(size_t spin_index) override {
        std::shared_lock lock(_mutex);
        if (spin_index >= this->_spin2cell->size()) {
            throw std::out_of_range(std::format("Spin index out of range (<{}>)", this->_spin2cell->size()));
        }
        return this->_spin2cell->at(spin_index);
    };

    /**
     * @brief Create macrocells if not already present, or recreate if requested.
     *
     * @param recreate Force regeneration of macrocell structure (default: false).
     *
     * @returns void - mutates internal state if needed.
     */
    PYTHON_API virtual void createMacrocellsIfNotCreated(bool recreate = false) override {
        std::unique_lock lock(_mutex);
        if (this->_macrocells->empty() || recreate) {
            this->createMacrocellsImpl(); // создать макроячейки, поделив спины
        }
    };

    /**
     * @brief Share external macrocell caches with another geometry.
     *
     * Used for cloning or caching scenarios.
     *
     * @param macrocell_index_cache Shared neighbor cache.
     * @param macrocells            Shared macrocell list.
     * @param spin2cell             Shared spin-to-cell mapping.
     *
     * @returns void - mutates internal pointers.
     */
    void shareMacrocellCache(
        std::shared_ptr<MacrocellIndexCache> macrocell_index_cache,
        std::shared_ptr<std::vector<NamespaceMacroCell>> macrocells,
        std::shared_ptr<std::vector<size_t>> spin2cell
    ) {
        std::unique_lock lock(_mutex);
        this->_macrocell_index_cache = macrocell_index_cache;
        this->_macrocells = macrocells;
        this->_spin2cell = spin2cell;
    };

    /**
     * @brief Get linear macrocell size (length of side) in meters.
     *
     * @returns Edge length of macrocell in SI units [m].
     */
    virtual double getMacrocellSize() const override { return this->macrocell_size; };

    /**
     * @brief Get macrocell volume in cubic meters.
     *
     * @returns Volume of a single macrocell in SI units [m³].
     */
    virtual double getMacrocellVolume() const override { return std::pow(this->macrocell_size, 3); };

    /**
     * @brief Retrieve neighboring macrocell moments within a cutoff radius.
     *
     * Performs spatial filtering of macrocells and returns their average moments.
     * Used in demagnetization/dipolar approximations for example.
     *
     * @param index          Spin index serving as query center.
     * @param cutoff_radius  Search radius [m].
     *
     * @returns Container of shared pointers to neighboring macrocell moments.
     */
    virtual MomentsContainer<CartesianCoordSystem> getMomentsFromMacrocells(
        size_t index, double cutoff_radius
    ) override {
        //
        size_t my_cell_index = this->getMacrocellIndexBySpin(index);

        std::shared_lock read_cache_lock(_mutex);
        // проверка на нуль
        if (cutoff_radius <= 1e-15)
            return {}; // типичные размеры: 1e-10

        // ---- проверяем кэш ----
        auto cache_key = std::make_pair(index, cutoff_radius);
        if (this->_macrocell_index_cache->has(cache_key)) {
            return this->getFromCacheImpl(cache_key);
        }

        read_cache_lock.unlock(); // разблокируем мьютекс, чтобы не блокировать другие потоки

        // и заблокируем его на запись (другим запрещаем лезть сюда)
        std::unique_lock write_lock(_mutex);

        // ---- проверяем кэш (могли записать за время) ----
        if (this->_macrocell_index_cache->has(cache_key)) {
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
        for (size_t i = 0; i < _macrocells->size(); ++i) {
            if (my_cell_index == i)
                continue; // пропускаем свою ячейку
            const double distance =
                target_coords.getDistanceFrom(_macrocells->at(i).avg_moment->getCoordinates());
            if (distance <= cutoff_radius) {
                neighbors.push_back(i);
            }
        }

        // ---- сохраняем в кэш ----
        this->_macrocell_index_cache->update(cache_key, neighbors);

        // и возвращаем список моментов из индексов макроячеек в радиусе обрезки
        MomentsContainer<CartesianCoordSystem> moments_from_macrocells(neighbors.size(), nullptr);
        for (size_t macro_index = 0; macro_index < neighbors.size(); ++macro_index) {
            moments_from_macrocells[macro_index] = this->_macrocells->at(neighbors[macro_index]).avg_moment;
        }
        return moments_from_macrocells;
    };

    /**
     * @brief Recalculate average values for each macrocell.
     *
     * Updates properties (average position, direction, material, amplitude, etc.) for each cell.
     *
     * @returns void - mutates macrocell average moment state.
     */
    PYTHON_API virtual void updateMacrocells() override {
        SCOPED_LOG_TIMER_DEBUG("│  ├─ Updating macrocells");
        // блокировка мьютекса
        std::unique_lock lock(_mutex);

        if (!this->_moments || this->_moments->empty()) {
            return;
        }

        // Итерируем по существующим ячейкам
        for (auto &cell : *this->_macrocells) {
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

            // "общая" намагниченность насыщения (для определения магнитного веса макроячейки)
            Eigen::Vector3d total_moment_vector(0, 0, 0);

            // обход по всем моментам в ячейке
            for (size_t spin_idx : cell.moment_indices) {
                moment = (*this->_moments)[spin_idx];
                avg_coordinates_vector += moment->getCoordinates().asVector();
                avg_direction_vector += moment->getDirection().asVector();
                if ((material = moment->getMaterialSharedPtr())) {
                    material_counts[material]++;
                    total_moment_vector += moment->getDirection().asVector() *
                                           material->atomic_magnetic_saturation_magnetization;
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

            // наконец, обновляем момент в ячейке
            cell.avg_moment =
                std::make_shared<Moment>(avg_coordinates_vector, avg_direction_vector, predominant_material);
            // общий вес намагниченности в диполь-дипольном взаимодействии
            cell.avg_moment->amplitude =
                total_moment_vector.norm() / predominant_material->atomic_magnetic_saturation_magnetization;
        }
    };

    /**
     * @brief Clear macrocell state.
     *
     * @returns void - mutates internal macrocell containers.
     */
    PYTHON_API virtual void clearMacrocells() override {
        std::unique_lock lock(_mutex);
        this->clearMacrocellsImpl();
    };
};

/**
 * @class Geometry
 * @brief Concrete geometry class implemented in the Cartesian coordinate system.
 *
 * Provides thread-safe access to a collection of magnetic moments,
 *   their topological layout, neighbor relations, macrocell grouping, etc.
 *
 * @details Implements the `AbstractGeometry` interface and macrocell handling via `MacrocellManager`.
 *
 * Supports construction from raw moment containers or structured numpy arrays
 *   including material assignments.
 *
 * @details Thread-safety: all methods that mutate internal state are protected by a shared mutex.
 *
 * Geometry instances own and expose:
 *   - All magnetic moments in the sample.
 *   - Spatial layout and neighbor information.
 *   - Macrocell-based spatial grouping and aggregation.
 *   - Read-only or cloned access to internal state.
 *
 * @note Units are SI unless otherwise stated.
 */
class Geometry : public AbstractGeometry, public MacrocellManager {
  protected:
    /**
     * @brief Container of magnetic moments in the geometry.
     *
     * Each moment is stored as a shared pointer. All moment-related queries
     *   (coordinates, directions, materials) resolve through this container.
     */
    MomentsContainer<NamespaceCoordSystem> _moments;

    /**
     * @brief Cache for neighbor indices per cutoff radius.
     *
     * Used to accelerate repeated spatial queries.
     * Keyed by (moment index, cutoff radius).
     */
    std::shared_ptr<MomentIndexCache> _moment_index_cache;

    /**
     * @brief Shared mutex guarding access to internal data structures.
     *
     * Allows concurrent reads and safe mutation.
     */
    mutable std::shared_mutex _mutex;

  public:
    /**
     * @brief Construct geometry from an existing moment container (Lvalue).
     *
     * @param moments         Container of shared pointers to moment objects.
     * @param macrocell_size  Macrocell discretization size [m].
     *
     * @throws std::invalid_argument
     */
    PYTHON_API Geometry(const MomentsContainer<NamespaceCoordSystem> &moments, double macrocell_size = 1e-9)
        : MacrocellManager(nullptr, macrocell_size),
          _moments(moments),
          _moment_index_cache(std::make_shared<MomentIndexCache>()) {
        MacrocellManager::_moments = &this->_moments; // заполнить указатель на массив моментов
    };

    /**
     * @brief Construct geometry from an existing moment container (Rvalue).
     *
     * @param moments         Container of shared pointers to moment objects.
     * @param macrocell_size  Macrocell discretization size [m].
     *
     * @throws std::invalid_argument
     */
    PYTHON_API Geometry(MomentsContainer<NamespaceCoordSystem> &&moments, double macrocell_size = 1e-9)
        : MacrocellManager(nullptr, macrocell_size),
          _moments(std::move(moments)),
          _moment_index_cache(std::make_shared<MomentIndexCache>()) {
        MacrocellManager::_moments = &this->_moments; // заполнить указатель на массив моментов
    };

    /**
     * @brief Construct geometry from a numpy-style matrix.
     *
     * Matrix must contain 7 columns per row: [x, y, z, sx, sy, sz, material_id].
     * Material registry is queried per row to construct correct material binding.
     *
     * @param moments           Matrix of shape (N, 7) representing the geometry.
     * @param material_registry Registry providing materials by integer ID.
     * @param macrocell_size    Macrocell discretization size [m].
     *
     * @throws std::invalid_argument
     */
    PYTHON_API Geometry(
        const Eigen::MatrixXd &moments, MaterialRegistry &material_registry, double macrocell_size = 1e-9
    )
        : MacrocellManager(nullptr, macrocell_size),
          _moment_index_cache(std::make_shared<MomentIndexCache>()) {
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

    /**
     * @brief Get number of magnetic moments in the geometry.
     *
     * @returns Number of moments (dimension of the system).
     */
    PYTHON_API virtual size_t size() const override {
        std::shared_lock lock(_mutex);
        return this->_moments.size();
    }

    /**
     * @brief Get a copy of all moments in the geometry.
     *
     * Returns a container of shared pointers to internal moment objects.
     * Can be used for analysis or iteration.
     *
     * @returns Container of all moments in the system.
     */
    PYTHON_API virtual MomentsContainer<NamespaceCoordSystem> getMoments() const override {
        std::shared_lock lock(_mutex);
        return this->_moments;
    };

    /**
     * @brief Optional geometry preparation step.
     *
     * This may construct internal spatial structures, create macrocells,
     *   or cache neighbor relations. May mutate internal state.
     * @note Will be replaced by proper constructor-based initialization in future.
     *
     * @returns void - mutates geometry and macrocell state.
     */
    virtual void prepare() override {
        std::unique_lock lock(_mutex); // мутирует _moments по указателю из менеджера
        this->createMacrocellsIfNotCreated(true); // создать макроячейки, поделив спины
        this->updateMacrocells();                 // посчитать их в первый раз
    };

    /**
     * @brief Access a mutable reference to a moment by index.
     *
     * @param index Index of the moment in the geometry container.
     *
     * @returns Reference to mutable moment.
     */
    PYTHON_API virtual Moment &operator[](size_t index) override {
        std::shared_lock lock(_mutex); // TODO может быть медленно!!!
        return *this->_moments[index];
    };

    /**
     * @brief Extract moments from specified indices.
     *
     * @param indexes Vector of moment indices to extract.
     *
     * @returns Container of corresponding moments.
     */
    virtual MomentsContainer<NamespaceCoordSystem> getFromIndexes(const std::vector<size_t> &indexes
    ) const override {
        std::shared_lock lock(_mutex);
        MomentsContainer<NamespaceCoordSystem> result(indexes.size(), nullptr);
        for (size_t i = 0; i < indexes.size(); ++i) {
            result[i] = this->_moments[indexes[i]];
        }
        return result;
    };

    /**
     * @brief Get neighboring moment indices for a given moment.
     *
     * Performs spatial filtering using Euclidean distance.
     *
     * @param index         Index of the center moment.
     * @param cutoff_radius Cutoff radius for neighbor search [m].
     *
     * @returns Vector of neighbor indices.
     */
    PYTHON_API virtual std::vector<size_t> getNeighbors(size_t index, double cutoff_radius) override {
        std::shared_lock lock(_mutex);
        // проверка на нуль
        if (cutoff_radius <= 1e-15)
            return {}; // типичные размеры: 1e-10

        // ---- проверяем кэш ----
        auto cache_key = std::make_pair(index, cutoff_radius);
        if (this->_moment_index_cache->has(cache_key)) {
            return this->_moment_index_cache->get(cache_key);
        }

        // ---- иначе честно считаем ----
        std::vector<size_t> neighbors;
        const auto &target_coords = this->_moments[index]->getCoordinates();

        for (size_t i = 0; i < this->_moments.size(); ++i) {
            if (i == index)
                continue; // себя не считаем
            const double distance = target_coords.getDistanceFrom(_moments[i]->getCoordinates());
            if (distance <= cutoff_radius) {
                neighbors.push_back(i);
            }
        }

        // ---- сохраняем в кэш ----
        this->_moment_index_cache->update(cache_key, neighbors);

        // возвращаем индексы соседей
        return neighbors;
    };

    // ---------------------------------------------------------------------
    //  Iterators
    // ---------------------------------------------------------------------

    /**
     * @brief Type alias for iterator over internal moments.
     */
    using iterator = MomentsContainer<NamespaceCoordSystem>::iterator;

    /**
     * @brief Iterator to the beginning of the moment container.
     *
     * @returns Iterator to first moment.
     */
    virtual iterator begin() override {
        std::shared_lock lock(_mutex);
        return this->_moments.begin();
    }

    /**
     * @brief Iterator to the end of the moment container.
     *
     * @returns Iterator past the last moment.
     */
    virtual iterator end() override {
        std::shared_lock lock(_mutex);
        return this->_moments.end();
    }

    /**
     * @brief Share spatial neighbor cache with another geometry.
     *
     * @param cache Shared pointer to neighbor index cache.
     *
     * @returns void - mutates internal state.
     */
    void shareMomentsIndexCache(const std::shared_ptr<MomentIndexCache> &cache) {
        std::unique_lock lock(this->_mutex);
        this->_moment_index_cache = cache;
    }

    /**
     * @brief Clone this geometry and return a new instance.
     *
     * Can optionally share internal caches and trigger preparation steps.
     *
     * @param share_cache Whether to reuse internal caches from the original.
     * @param prepare     Whether to invoke `prepare()` after cloning.
     *
     * @returns Unique pointer to new geometry instance.
     */
    virtual std::unique_ptr<IGeometry<NamespaceCoordSystem>> clone(
        bool share_cache = false, bool prepare = false
    ) const override {
        std::shared_lock lock(_mutex);
        MomentsContainer<NamespaceCoordSystem> cloned_moments;
        cloned_moments.reserve(this->_moments.size());
        for (const auto &moment : this->_moments) {
            cloned_moments.push_back(moment->clone());
        }
        auto cloned = std::make_unique<Geometry>(cloned_moments, this->macrocell_size);
        if (share_cache) {
            cloned->shareMomentsIndexCache(this->_moment_index_cache);
            cloned->shareMacrocellCache(this->_macrocell_index_cache, this->_macrocells, this->_spin2cell);
        }
        if (prepare) {
            cloned->prepare();
        }
        return cloned;
    }

    // ---------------------------------------------------------------------
    //  Representations
    // ---------------------------------------------------------------------

    /**
     * @brief Export geometry as NumPy-like matrix.
     *
     *  Each row contains: [x, y, z, sx, sy, sz, material_id].
     *
     * @returns Matrix of shape (N, 7).
     */
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

    /**
     * @brief Get human-readable string representation.
     *
     * Typically used for debugging or logging.
     *
     * @returns Stringified geometry contents.
     */
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

// ===========================================================================
//  Python bindings
// ===========================================================================

#define MACROCELL_TEMPLATE_BINDINGS(cls)                                                                     \
    .def_readonly(                                                                                           \
        "moment_indices",                                                                                    \
        &cls::moment_indices,                                                                                \
        py::doc(                                                                                             \
            "@brief Indices of moments contained within this macrocell.\n"                                   \
            "\n"                                                                                             \
            "Represents the positions (indices) in the global geometry container that belong to this cell."  \
        )                                                                                                    \
    )                                                                                                        \
        .def_property_readonly(                                                                              \
            "avg_moment",                                                                                    \
            [](const cls &self) { return self.avg_moment.get(); },                                           \
            py::return_value_policy::reference_internal,                                                     \
            py::doc("@brief Average moment of the macrocell.\n"                                              \
                    "\n"                                                                                     \
                    "Computed from all internal moments using vector averaging of position and direction.\n" \
                    "@details Used in demagnetization, dipole-dipole interactions, etc.")                    \
        )

#define MACROCELLMANAGER_TEMPLATE_BINDINGS(cls)                                                              \
    .def(                                                                                                    \
        "get_macrocells",                                                                                    \
        &cls::getMacrocells,                                                                                 \
        py::doc("@brief Access the full macrocell container.\n"                                              \
                "\n"                                                                                         \
                "@returns Reference to internal vector of macrocells.")                                      \
    )                                                                                                        \
        .def(                                                                                                \
            "get_macrocell_index_by_spin",                                                                   \
            &cls::getMacrocellIndexBySpin,                                                                   \
            py::doc("@brief Get the macrocell index for a given spin index.\n"                               \
                    "\n"                                                                                     \
                    "Maps a spin index to its owning macrocell.\n"                                           \
                    "@param spin_index Index of spin in global geometry.\n"                                  \
                    "\n"                                                                                     \
                    "@returns Index of macrocell containing this spin.\n")                                   \
        )                                                                                                    \
        .def(                                                                                                \
            "update_macrocells",                                                                             \
            &cls::updateMacrocells,                                                                          \
            py::doc("@brief Recalculate average values for each macrocell.\n"                                \
                    "\n"                                                                                     \
                    "Updates properties (average position, direction, material, amplitude, etc.) for each "  \
                    "cell.\n"                                                                                \
                    "\n"                                                                                     \
                    "@returns void - mutates macrocell average moment state.\n")                             \
        )                                                                                                    \
        .def(                                                                                                \
            "create_macrocells_if_not_created",                                                              \
            &cls::createMacrocellsIfNotCreated,                                                              \
            py::arg("recreate") = false,                                                                     \
            py::doc("@brief Create macrocells if not already present, or recreate if requested.\n"           \
                    "\n"                                                                                     \
                    "@param recreate Force regeneration of macrocell structure (default: false).\n"          \
                    "\n"                                                                                     \
                    "@returns void - mutates internal state if needed.\n")                                   \
        )                                                                                                    \
        .def(                                                                                                \
            "clear_macrocells",                                                                              \
            &cls::clearMacrocells,                                                                           \
            py::doc("@brief Clear macrocell state.\n"                                                        \
                    "\n"                                                                                     \
                    "@returns void - mutates internal macrocell containers.\n")                              \
        )

#define GEOMETRY_TEMPLATE_BINDINGS(cls)                                                                      \
    .def(                                                                                                    \
        "__str__",                                                                                           \
        &cls::__str__,                                                                                       \
        py::doc("@brief Get human-readable string representation.\n"                                         \
                "\n"                                                                                         \
                "Typically used for debugging or logging.\n"                                                 \
                "@returns Stringified geometry contents.")                                                   \
    )                                                                                                        \
        .def(                                                                                                \
            "__len__",                                                                                       \
            &cls::size,                                                                                      \
            py::doc("@brief Get number of magnetic moments in the geometry.\n"                               \
                    "\n"                                                                                     \
                    "@returns Number of moments (dimension of the system).")                                 \
        )                                                                                                    \
        .def(                                                                                                \
            "get_neighbors",                                                                                 \
            &cls::getNeighbors,                                                                              \
            py::arg("index"),                                                                                \
            py::arg("cutoff_radius"),                                                                        \
            py::doc("@brief Get neighboring moment indices for a given moment.\n"                            \
                    "\n"                                                                                     \
                    "Performs spatial filtering using Euclidean distance.\n\n"                               \
                    "@param index         Index of the center moment.\n"                                     \
                    "@param cutoff_radius Cutoff radius for neighbor search [m].\n"                          \
                    "\n"                                                                                     \
                    "@returns Vector of neighbor indices.")                                                  \
        )                                                                                                    \
        .def(                                                                                                \
            "as_numpy",                                                                                      \
            &cls::asNumpy,                                                                                   \
            py::doc("@brief Export geometry as NumPy-like matrix.\n"                                         \
                    "\n"                                                                                     \
                    "Each row contains:\n"                                                                   \
                    "  coords, direction, material, etc.\n"                                                  \
                    "\n"                                                                                     \
                    "@returns Matrix of shape (N, xxx).\n")                                                  \
        )                                                                                                    \
        .def(                                                                                                \
            "__getitem__",                                                                                   \
            &cls::operator[],                                                                                \
            py::arg("index"),                                                                                \
            py::return_value_policy::reference,                                                              \
            py::doc("@brief Access a mutable reference to a moment by index.\n"                              \
                    "\n"                                                                                     \
                    "@param index Index of the moment in the geometry container.\n"                          \
                    "\n"                                                                                     \
                    "@returns Reference to mutable moment.\n")                                               \
        )                                                                                                    \
        .def(                                                                                                \
            "get_moments",                                                                                   \
            &cls::getMoments,                                                                                \
            py::doc("@brief Get a copy of all moments in the geometry.\n"                                    \
                    "\n"                                                                                     \
                    "Returns a container of shared pointers to internal moment objects.\n"                   \
                    "Can be used for analysis or iteration.\n"                                               \
                    "\n"                                                                                     \
                    "@returns Container of all moments in the system.\n")                                    \
        )

/**
 * @brief Bind the geometry utilities to a Python sub‑module.
 *
 * @param module Parent PyBind11 module (usually the core extension module).
 * @returns void – extends the parent module.
 */
inline void pyBindGeometries(py::module_ &module) {
    using namespace spindynapy;

    // -------- | GEOMETRIES | --------
    py::module_ geometries_module = module.def_submodule("geometries");

    geometries_module.doc() =
        "@brief  Geometry interfaces and macrocell management for spin dynamics simulations.\n"
        "\n"
        "This header defines the abstract interfaces and concrete implementations\n"
        "  for spatial geometry and macrocell management used in spin dynamics simulations.\n"
        "\n"
        "Functional overview:\n"
        "\n"
        "- A **Geometry** manages the spatial layout of magnetic moments in a sample.\n"
        "  It provides access to moment positions, directions, and materials.\n"
        "  It also allows neighborhood queries and macrocell-based approximations.\n"
        "\n"
        "- A **MacrocellManager** subdivides the geometry into spatial regions (\" macrocells \"),\n"
        "  which aggregate moments for use in demagnetization and dipolar field calculations.\n"
        "\n"
        "- Thread safety:\n"
        "  - All geometry implementations (e.g., `Geometry` in Cartesian coordinates)\n"
        "    are thread-safe via shared mutexes.\n"
        "  - Neighbor and macrocell lookups are cached using thread-safe caches.\n"
        "\n"
        "Interfaces:\n"
        "- `IGeometry<CoordSystem>`         – interface for geometry access and manipulation,\n"
        "- `IMacrocellManager<CoordSystem>`– interface for macrocell-based spatial aggregation,\n"
        "- `MacroCell<CoordSystem>`        – container for grouped moment indices and average moment,\n"
        "- `MomentsContainer<CoordSystem>` – alias for list of shared pointers to moments.\n"
        "\n"
        "Concrete (Cartesian):\n"
        "- `Geometry`                       – thread-safe geometry class in Cartesian system,\n"
        "- `MacrocellManager`              – thread-safe manager for Cartesian macrocells,\n"
        "- `MomentIndexCache`              – cache for moment neighbors by radius,\n"
        "- `MacrocellIndexCache`           – cache for macrocell neighbors by radius.\n"
        "\n"
        "Python bindings:\n"
        "- Macros:\n"
        "  - `MACROCELL_TEMPLATE_BINDINGS`\n"
        "  - `MACROCELLMANAGER_TEMPLATE_BINDINGS`\n"
        "  - `GEOMETRY_TEMPLATE_BINDINGS`\n"
        "- Binder function:\n"
        "  - `pyBindGeometries(py::module_ &)`\n"
        "\n"
        "@note All units are SI unless explicitly stated. Coordinates are in meters.\n"
        "      Directions are dimensionless unit vectors. Fields are in [T], energies in [J].\n"
        "\n"
        "@note This header only covers Cartesian coordinates; other systems may extend it in the future.\n";

    // -------- | CARTESIAN GEOMETRIES | --------
    py::module_ cartesian = geometries_module.def_submodule("cartesian");

    cartesian.doc() = geometries_module.doc();
    {
        using cartesian::AbstractGeometry;
        using cartesian::AbstractMacrocellManager;
        using cartesian::Geometry;
        using cartesian::MacrocellManager;
        using cartesian::Moment;
        using cartesian::NamespaceMacroCell;

        py::class_<NamespaceMacroCell, std::shared_ptr<NamespaceMacroCell>>(cartesian, "MacroCell")
            MACROCELL_TEMPLATE_BINDINGS(NamespaceMacroCell)
                .doc() =
            "@struct MacroCell\n"
            "@brief Macrocell container structure for grouping magnetic moments.\n"
            "\n"
            "Represents a logical subdivision of the geometry into macrocells. Each macrocell aggregates a "
            "subset\n"
            "  of magnetic moments and stores their average moment for use in approximated long-range "
            "interactions\n"
            "  such as demagnetization or dipole-dipole coupling.\n"
            "\n"
            "The macrocell contains:\n"
            "  - A list of indices of moments it owns.\n"
            "  - A shared pointer to the average magnetic moment computed from all internal moments.\n"
            "\n"
            "@note The average moment is updated via `updateMacrocells()` and is used in field "
            "approximations.\n";

        py::class_<AbstractMacrocellManager, std::shared_ptr<AbstractMacrocellManager>>(
            cartesian, "AbstractMacrocellManager"
        ) MACROCELLMANAGER_TEMPLATE_BINDINGS(AbstractMacrocellManager)
            .doc() =
            "@class AbstractMacrocellManager\n"
            "@brief Abstract interface for macrocell management in a geometry.\n"
            "\n"
            "This interface defines operations for managing macrocells in a geometry (moments lists)\n"
            "  for a given coordinate system. Macrocells are logical groups of spins\n"
            "  used to accelerate long-range interactions via spatial aggregation.\n"
            "\n"
            "Responsibilities include:\n"
            " - Macrocell creation and re-creation.\n"
            " - Updating average spin states within cells.\n"
            " - Mapping between spins and their macrocells.\n"
            " - Retrieving moments within nearby macrocells\n"
            "etc.\n"
            "\n"
            "@note This interface is extended by concrete geometry implementations and provides\n"
            "      internal mechanisms to interact with the macrocell cache.\n"
            "\n"
            "@see MacroCell\n";

        py::class_<AbstractGeometry, AbstractMacrocellManager, std::shared_ptr<AbstractGeometry>>(
            cartesian, "AbstractGeometry"
        ) GEOMETRY_TEMPLATE_BINDINGS(AbstractGeometry)
            .doc() = "@class AbstractGeometry\n"
                     "@brief Abstract interface for geometry in a cartesian coordinate system.\n"
                     "\n"
                     "This interface defines operations for accessing and manipulating\n"
                     "  the physical layout of magnetic moments in a system. It also extends\n"
                     "  macrocell-related functionality from `IMacrocellManager`.\n"
                     "\n"
                     "Geometry instances own and expose:\n"
                     "  - All magnetic moments in the sample.\n"
                     "  - Spatial layout and neighbor information.\n"
                     "  - Macrocell-based spatial grouping and aggregation.\n"
                     "  - Read-only or cloned access to internal state.\n"
                     "\n"
                     "@note All units are SI.\n";

        py::class_<MacrocellManager, AbstractMacrocellManager, std::shared_ptr<MacrocellManager>>(
            cartesian, "MacrocellManager"
        )
            .def_readonly(
                "macrocell_size",
                &MacrocellManager::macrocell_size,
                py::doc("@brief Size of each macrocell in all directions [m].\n"
                        "\n"
                        "@note Macrocell is a cubic structure, so this defines the edge length.\n"
                        "Defines the granularity of spatial decomposition.\n")
            )
            .doc() =
            "@class MacrocellManager\n"
            "@brief Thread-safe manager for spatial macrocell decomposition in Cartesian coordinate "
            "systems.\n"
            "\n"
            "Provides creation, caching, and access to macrocells, which group atomic spins into\n"
            "  spatial regions used for approximated long-range interactions such as demagnetization.\n"
            "\n"
            "Internally uses a regular grid aligned with the geometry bounding box, with\n"
            "  fixed macrocell size. Supports cache-aware neighbor lookup and shared cache usage.\n"
            "\n"
            "Responsibilities include:\n"
            " - Macrocell creation and re-creation.\n"
            " - Updating average spin states within cells.\n"
            " - Mapping between spins and their macrocells.\n"
            " - Retrieving moments within nearby macrocells\n"
            "etc.\n"
            "\n"
            "@note This class assumes a fixed-size cubic macrocell and uses a thread-safe shared cache\n"
            "  for neighbor lookups in the cutoff radius.\n";

        py::class_<Geometry, AbstractGeometry, MacrocellManager, std::shared_ptr<Geometry>>(
            cartesian, "Geometry"
        )
            .def(
                py::init<const MomentsContainer<CartesianCoordSystem> &, double>(),
                py::arg("moments"),
                py::arg("macrocell_size") = 1e-9
            )
            .def(
                py::init<const Eigen::MatrixXd &, MaterialRegistry &, double>(),
                py::arg("moments"),
                py::arg("material_registry"),
                py::arg("macrocell_size") = 1e-9
            )
            .doc() = "@class Geometry\n"
                     "@brief Concrete geometry class implemented in the Cartesian coordinate system.\n"
                     "\n"
                     "Provides thread-safe access to a collection of magnetic moments,\n"
                     "  their topological layout, neighbor relations, macrocell grouping, etc.\n"
                     "\n"
                     "@details Implements the `AbstractGeometry` interface and macrocell handling via "
                     "`MacrocellManager`.\n"
                     "\n"
                     "Supports construction from raw moment containers or structured numpy arrays\n"
                     "  including material assignments.\n"
                     "\n"
                     "@details Thread-safety: all methods that mutate internal state are protected by a "
                     "shared mutex.\n"
                     "\n"
                     "Geometry instances own and expose:\n"
                     "  - All magnetic moments in the sample.\n"
                     "  - Spatial layout and neighbor information.\n"
                     "  - Macrocell-based spatial grouping and aggregation.\n"
                     "  - Read-only or cloned access to internal state.\n"
                     "\n"
                     "@note Units are SI unless otherwise stated.\n";
    }
}

#endif // ! __GEOMETRIES_BASE_HPP__
