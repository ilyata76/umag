#ifndef __SPINDYNAPY_LOGGER_HPP__
#define __SPINDYNAPY_LOGGER_HPP__

/**
 * @file   logger.hpp
 * @brief  Diagnostic logging utility for SpinDynaPy (thread‑safe singleton and scoped timer).
 *
 * This header implements a lightweight, header‑only facility that lets the C++ core emit
 *   timestamped UTF‑8 logging. Messages are buffered in‑memory and written
 *   to a user‑supplied `std::ostream` only when `flush()` is requested,
 *   thereby amortising I/O cost across many log events.
 *
 * Exposed entities
 * - `spindynapy::Logger`      – global singleton for buffered logging.
 * - `spindynapy::ScopedTimer` – RAII helper that logs the execution time of a scope.
 * - Convenience macros: `LOG_MSG`, `LOG_MSG_PRINT`, `SCOPED_LOG_TIMER`, `SCOPED_LOG_TIMER_PRINT`,
 *      `SCOPED_LOG_TIMER_DEBUG`
 * - Function `pyBindLogger()` that exports the facility to Python as sub‑module `core.logger`.

 * @note The logger is intended for debugging and profiling. It is **not** a persistent audit trail.
 *
 * @copyright 2025 SpinDynaPy
 */

#include "constants.hpp"

#include <chrono>
#include <deque>
#include <iostream>
#include <mutex>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <string>

// ===========================================================================
//  Internal helper
// ===========================================================================

/**
 * @brief Obtain the current local time as a printable timestamp.
 *
 * The function is header‑only and `noexcept`; the formatting buffer is fixed‑size so no dynamic
 *   allocation occurs. Platform‑specific thread‑safe variants (`localtime_r`, `localtime_s`) are
 *   selected at compile time.
 *
 * @returns A string formatted as `YYYY‑MM‑DD HH:MM:SS`.
 */
inline std::string iso_timestamp_now() noexcept {
    using namespace std::chrono;
    const auto tp = system_clock::now();
    std::time_t t = system_clock::to_time_t(tp);
    std::tm tm;
#if defined(_MSC_VER)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    char buf[20]; // 19 chars + null‑terminator
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &tm);
    return buf;
}

namespace spindynapy {

// ===========================================================================
//  LOGGER
// ===========================================================================

/**
 * @class   Logger
 * @brief   Thread‑safe singleton that buffers log messages and flushes them on demand.
 *
 * The logger prepends every message with a human‑readable timestamp.
 * Messages are stored in a `std::deque` to minimise allocation churn and are written
 *   to the configured output stream **only** when `flush()` is invoked,
 *   thereby amortising I/O cost across many log events.
 *
 * Guarantees
 * - Singleton access via `Logger::instance()`.
 * - All mutating operations (`add`, `flush`, `setStream`, `resetStream`) are mutex‑protected.
 * - Default output stream is `std::cout`; users may redirect to any `std::ostream` that outlives the logger
 *
 * @note The singleton is never destroyed; relying modules may therefore safely log from static
 *         destructors.
 */
class PYTHON_API Logger {
  private:
    /** @brief Private constructor; use `instance()` instead. */
    Logger() noexcept : _out(&std::cout) {}

    /** Output destination – defaults to `std::cout`. (not owned) */
    std::ostream *_out;

    /** Buffer of fully‑formatted log lines. */
    std::deque<std::string> _buffer;

    /** Mutex protecting both the buffer and the output pointer. */
    std::mutex _mutex;

  public:
    // Deleted copy / move to enforce singleton semantics.
    Logger(const Logger &) = delete;
    // Deleted copy / move to enforce singleton semantics.
    Logger &operator=(const Logger &) = delete;
    // Deleted copy / move to enforce singleton semantics.
    Logger(Logger &&) = delete;
    // Deleted copy / move to enforce singleton semantics.
    Logger &operator=(Logger &&) = delete;

    /**
     * @brief Retrieve the global logger instance.
     *
     * First call lazily initialises the static instance.
     * Subsequent calls are inexpensive reference returns.
     *
     * @returns Reference to the singleton logger.
     */
    PYTHON_API static Logger &instance() noexcept {
        static Logger static_logger;
        return static_logger;
    }

    /**
     * @brief Route log output to a user‑provided stream.
     *
     * @note Thread‑safe
     *
     * @param os Destination output stream that must remain valid until replaced
     *              or `reset_to_stdout()` is invoked.
     * @returns void. Mutates internal output parameter by provided stream.
     */
    PYTHON_API void setStream(std::ostream &os) noexcept {
        std::lock_guard<std::mutex> lock(this->_mutex);
        this->_out = &os;
    }

    /**
     * @brief Convenience wrapper restoring logging to `stdout`.
     *
     * @details Reset to `std::cout`;
     * @note Thread‑safe
     *
     * @returns void. Mutates internal output parameter to default.
     */
    PYTHON_API void resetStream() noexcept { this->setStream(std::cout); }

    /**
     * @brief Append a message to the internal buffer.
     *
     * The message is timestamped immediately; no I/O occurs until `flush()` is called.
     * @note Thread‑safe
     *
     * @param msg UTF‑8 string to log.
     * @returns void – adds one element to the buffer.
     */
    PYTHON_API void add(const std::string &msg) noexcept {
        std::lock_guard<std::mutex> lock(this->_mutex);
        this->_buffer.emplace_back("[" + iso_timestamp_now() + "] " + msg);
    }

    /**
     * @brief Write all buffered messages to the configured stream and clear the buffer.
     *
     * @note Thread‑safe
     *
     * @returns void – empties the buffer and flushes the stream.
     */
    PYTHON_API void flush() {
        std::lock_guard<std::mutex> lock(this->_mutex);
        while (!this->_buffer.empty()) {
            (*this->_out) << this->_buffer.front() << '\n';
            this->_buffer.pop_front();
        }
        this->_out->flush();
    }
};

// ===========================================================================
//  ScopedTimer
// ===========================================================================

/**
 * @class   ScopedTimer
 * @brief   RAII helper that logs the duration of a scope in microseconds.
 *
 * The timer logs a `[START]` message upon construction and an `[END]` message upon destruction
 *   containing the elapsed time.
 * It is intended for quick instrumentation of code blocks without manual time bookkeeping.
 *
 * Usage example:
 * @code{.cpp} {
 *     {
 *          spindynapy::ScopedTimer timer("matrix‑multiplication");
 *          multiply(A, B, C);
 *     }
 * } @endcode
 *
 * @note Thread‑safe as long as the underlying `Logger` implementation is thread‑safe.
 */
class PYTHON_API ScopedTimer {
  private:
    /** Logical name of the timed scope. */
    std::string _name;

    /** Logger POINTER used for output (may be null). (not owned) */
    Logger *_logger;

    /** Whether to flush the logger immediately after each message. */
    bool _always_flush;

    /** Start time recorded at construction. */
    std::chrono::high_resolution_clock::time_point _start;

  public:
    /**
     * @brief Construct a timer and emit the start message to logger.
     *
     * @param name   Logical name of the timed scope (will appear in log).
     * @param logger Logger instance to use (defaults to global singleton).
     * @returns void – constructs and records start time.
     */
    PYTHON_API ScopedTimer(
        const std::string &name, bool always_flush = false, Logger *logger = &Logger::instance()
    )
        : _name(name), _logger(logger), _start(std::chrono::high_resolution_clock::now()) {
        if (this->_logger) {
            this->_logger->add("[START] " + _name);
            this->_always_flush = always_flush;
            if (this->_always_flush)
                this->_logger->flush(); // Flush immediately to ensure visibility
        }
    }

    /**
     * @brief Destruct the timer and emit to logger the end message with duration.
     *
     * @returns void – mutates logger state by adding an `[END]` entry.
     */
    PYTHON_API ~ScopedTimer() {
        if (this->_logger) {
            auto end = std::chrono::high_resolution_clock::now();
            auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - _start).count();
            this->_logger->add(
                "  [END] " + _name + " | duration: " + std::format(std::locale(""), "{:L}", dur) + " μs"
            );
            if (this->_always_flush)
                this->_logger->flush(); // Flush immediately to ensure visibility
        }
    }
};

// ===========================================================================
//  Convenience macros
// ===========================================================================

/** @brief Append MSG to global logger.
 *  @param MSG UTF‑8 text to log.
 *  @returns void – enqueues message. */
#define LOG_MSG(msg) ::spindynapy::Logger::instance().add(msg);

/** @brief Log MSG and immediately flush.
 *  @param MSG UTF‑8 text to log.
 *  @returns void – enqueues message and performs I/O. */
#define LOG_MSG_PRINT(msg)                                                                                   \
    ::spindynapy::Logger::instance().add(msg);                                                               \
    ::spindynapy::Logger::instance().flush();

/** @brief Create a scoped timer named NAME.
 *  @param NAME Logical name of the timed scope. */
#define SCOPED_LOG_TIMER(name)                                                                               \
    ::spindynapy::ScopedTimer _CONCAT(_scoped_timer_, __LINE__)(name, &Logger::instance());

#ifdef DEBUG
/** @brief Create a scoped timer named NAME.
 *  @param NAME Logical name of the timed scope. */
#define SCOPED_LOG_TIMER_DEBUG(name)                                                                         \
    ::spindynapy::ScopedTimer _CONCAT(_scoped_timer_, __LINE__)(name, &Logger::instance());
#else
/** @brief Create a scoped timer named NAME.
 *  @param NAME Logical name of the timed scope. */
#define SCOPED_LOG_TIMER_DEBUG(name)                                                                         \
    ::spindynapy::ScopedTimer _CONCAT(_scoped_timer_, __LINE__)(name, &Logger::instance());
#endif

/** @brief Create a scoped timer named NAME and flush when the scope ends.
 *  @param NAME Logical name of the timed scope. */
#define SCOPED_LOG_TIMER_PRINT(name)                                                                         \
    ::spindynapy::ScopedTimer _CONCAT(_scoped_timer_, __LINE__)(name, &Logger::instance());                  \
    ::spindynapy::Logger::instance().flush();

} // namespace spindynapy

// ===========================================================================
//  Python bindings
// ===========================================================================

/**
 * @brief Bind the logging utilities to a Python sub‑module.
 *
 * @param module Parent PyBind11 module (usually the core extension module).
 * @returns void – extends the parent module.
 */
inline void pyBindLogger(py::module_ &module) {
    using namespace spindynapy;

    // ---------------------------------------------------------------------
    //  Singleton specifics
    // ---------------------------------------------------------------------

    /**
     * @brief Adapter that exposes a Python file-like object as `std::ostream`.
     *
     * Holds a strong reference (`keep_alive`) so the Python object remains alive while C++
     *   writes through the STL interface (`os`).  Internally relies on
     *   `py::detail::pythonbuf` to bridge `ostream` and the Python `write()` protocol.
     *
     * @note Non-copyable, internal use only.
     */
    struct IOStreamAdapter {
        py::object keep_alive;                      ///< Python stream (ref-counted).
        std::unique_ptr<py::detail::pythonbuf> buf; ///< pybind11 buffer wrapper.
        std::ostream os;                            ///< Facade over `buf`.
        explicit IOStreamAdapter(const py::object &o) noexcept
            : keep_alive(o), buf(std::make_unique<py::detail::pythonbuf>(keep_alive)), os(buf.get()) {}

        /**
         * @brief Destructor releases the Python buffer and resets the stream.
         * @details If Python is not initialized, the buffer is leaked intentionally
         *          to avoid dangling references.
         */
        ~IOStreamAdapter() noexcept {
            os.rdbuf(nullptr);
            if (Py_IsInitialized()) {
                py::gil_scoped_acquire gil;
                buf.reset();
            } else {
                [[maybe_unused]] auto l =
                    buf.release(); // Release ownership if Python is not initialized
                                   // (leakage is acceptable here as this is a singleton).
            }
        }
    };

    /** @brief Non-owning pointer type for exposing the C++ singleton to Python. */
    using NoDelete = std::unique_ptr<Logger, py::nodelete>;

    /**
     * @brief Keeps the last Python stream alive while the logger uses it.
     *
     * When `Logger.set_stream()` is called from Python we wrap the object in an
     * `IOStreamAdapter` and store it here so it cannot be garbage-collected.
     */
    static std::shared_ptr<IOStreamAdapter> singleton_logger_holder;

    // ---------------------------------------------------------------------
    //  Create submodule
    // ---------------------------------------------------------------------

    py::module_ logger_module = module.def_submodule("logger");

    logger_module.doc() =
        ("@file   logger.hpp\n"
         "@brief  Diagnostic logging utility for SpinDynaPy (thread-safe singleton and scoped timer).\n"
         "\n"
         "This header implements a lightweight, header-only facility that lets the C++ core emit\n"
         "  timestamped UTF-8 logging. Messages are buffered in-memory and written\n"
         "  to a user-supplied `std::ostream` only when `flush()` is requested,\n"
         "  thereby amortising I/O cost across many log events.\n"
         "\n"
         "Exposed entities\n"
         "- `spindynapy::Logger`      – global singleton for buffered logging.\n"
         "- `spindynapy::ScopedTimer` – RAII helper that logs the execution time of a scope.\n"
         "- Convenience macros: `LOG_MSG`, `LOG_MSG_PRINT`, `SCOPED_LOG_TIMER`, `SCOPED_LOG_TIMER_PRINT`.\n"
         "- Function `pyBindLogger()` that exports the facility to Python as sub-module `core.logger`.\n"
         "\n"
         "@note The logger is intended for debugging and profiling. It is **not** a persistent audit trail.\n"
         "\n"
         "@copyright 2025 SpinDynaPy");

    // ---------------------------------------------------------------------
    //  Logger class binding
    // ---------------------------------------------------------------------

    py::class_<Logger, NoDelete>(logger_module, "Logger")
        .def(py::init([]() { return &Logger::instance(); }), py::return_value_policy::reference)
        .def_static(
            "instance",
            []() -> Logger & { return Logger::instance(); },
            py::return_value_policy::reference,
            py::doc("@brief Retrieve the global logger instance.\n"
                    "\n"
                    "First call lazily initialises the static instance.\n"
                    "Subsequent calls are inexpensive reference returns.\n"
                    "\n"
                    "@returns Reference to the singleton logger.")
        )
        .def(
            "add",
            &Logger::add,
            py::arg("msg"),
            py::doc("@brief Append a message to the internal buffer.\n"
                    "\n"
                    "The message is timestamped immediately; no I/O occurs until `flush()` is called.\n"
                    "\n"
                    "@note Thread-safe\n"
                    "\n"
                    "@param msg UTF-8 string to log.\n"
                    "@returns void – adds one element to the buffer.")
        )
        .def(
            "flush",
            &Logger::flush,
            py::doc("@brief Write all buffered messages to the configured stream and clear the buffer.\n"
                    "\n"
                    "@note Thread-safe\n"
                    "\n"
                    "@returns void – empties the buffer and flushes the stream.")
        )
        .def(
            "set_stream",
            [](Logger &self, py::object stream) {
                if (stream.is_none()) {
                    singleton_logger_holder.reset();
                    self.resetStream();
                    return;
                }
                singleton_logger_holder = std::make_shared<IOStreamAdapter>(stream);
                self.setStream(singleton_logger_holder->os);
            },
            py::arg("stream"),
            py::doc("@brief Route log output to a user-provided stream.\n"
                    "\n"
                    "@note Thread-safe\n"
                    "\n"
                    "@param stream Destination output stream that must remain valid until replaced\n"
                    "              or `reset_to_stdout()` is invoked. Use `None` to restore stdout.\n"
                    "@returns void. Mutates internal output parameter by provided stream.")
        )

        .def(
            "reset_stream",
            [](Logger &self) {
                singleton_logger_holder.reset();
                self.resetStream();
            },
            py::doc("@brief Convenience wrapper restoring logging to `stdout`.\n"
                    "\n"
                    "@details Reset to `std::cout`;\n"
                    "@note Thread-safe\n"
                    "\n"
                    "@returns void. Mutates internal output parameter to default.")
        )
        .doc() =
        ("@brief   Thread-safe singleton that buffers log messages and flushes them on demand.\n"
         "\n"
         "The logger prepends every message with a human-readable timestamp.\n"
         "Messages are stored in a `std::deque` to minimise allocation churn and are written\n"
         "  to the configured output stream **only** when `flush()` is invoked,\n"
         "  thereby amortising I/O cost across many log events.\n"
         "\n"
         "Guarantees\n"
         " - Singleton access via `Logger::instance()`.\n"
         " - All mutating operations (`add`, `flush`, `setStream`, `resetStream`) are mutex-protected.\n"
         " - Default output stream is `std::cout`; users may redirect to any `std::ostream` that outlives "
         "the "
         "logger\n"
         "\n"
         "@note The singleton is never destroyed; relying modules may therefore safely log from static\n"
         "        destructors.");

    // ---------------------------------------------------------------------
    //  ScopedTimer binding
    // ---------------------------------------------------------------------

    py::class_<spindynapy::ScopedTimer>(logger_module, "ScopedTimer")
        .def(
            py::init<const std::string &, bool, spindynapy::Logger *>(),
            py::arg("name"),
            py::arg("always_flush") = false,
            py::arg("logger") = &Logger::instance()
        )
        .def(
            "__enter__",
            [](spindynapy::ScopedTimer &t) { return &t; },
            py::doc("@brief Support for Python context-manager enter, starting timer\n\n@returns self.")
        )
        .def(
            "__exit__",
            [](spindynapy::ScopedTimer &, py::args) {},
            py::doc("@brief Support for Python context-manager exit; timer is closed automatically.")
        )
        .doc() =
        ("@class   ScopedTimer\n"
         "@brief   RAII helper that logs the duration of a scope in microseconds.\n"
         "\n"
         "The timer logs a `[START]` message upon construction and an `[END]` message upon destruction\n"
         "  containing the elapsed time.\n"
         "It is intended for quick instrumentation of code blocks without manual time bookkeeping.\n"
         "\n"
         "Usage example:\n"
         "```cpp\n"
         "{\n"
         "    spindynapy::ScopedTimer timer(\"matrix-multiplication\");\n"
         "    multiply(A, B, C);\n"
         "}\n"
         "```\n"
         "\n"
         "@note Thread-safe as long as the underlying `Logger` implementation is thread-safe.");
}

#endif // ! __SPINDYNAPY_LOGGER_HPP__
