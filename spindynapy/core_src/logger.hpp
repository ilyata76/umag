#ifndef __LOGGER_HPP__
#define __LOGGER_HPP__

/**
 * Логирование сообщений в stdout или в любой другой поток вывода.
 *   Для отладки и мониторинга выполнения программы, длительности исполнения блоков кода.
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

// Строка с текущими датой‑времем «YYYY-MM-DD HH:MM:SS»
inline std::string now_timestamp_str() {
    using namespace std::chrono;
    auto tp = system_clock::now();
    std::time_t t = system_clock::to_time_t(tp);
    std::tm tm;
#if defined(_MSC_VER)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    char buf[20]; // «YYYY-MM-DD HH:MM:SS» = 19 символов + 0
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &tm);
    return buf;
}

namespace spindynapy {

/**
 * Класс для логирования сообщений в stdout или в любой другой поток вывода.
 * !!! СИНГЛТОН !!! (потокобезопасный)
 *
 * Используется для отладки и мониторинга выполнения программы.
 * Вызывается через Logger::instance()
 */
class PYTHON_API Logger {
  private:
    // приватный конструктор синглтона
    Logger() : _out(&std::cout) {}

    // поток вывода сообщений
    std::ostream *_out;
    // хранимый буфер в виде сообщений
    std::deque<std::string> _buffer;
    // мьютекс для потокобезопасности
    std::mutex _mutex;

  public:
    // Получить объект логгера (синглтон)
    PYTHON_API static Logger &instance() {
        static Logger static_logger;
        return static_logger;
    }

    // Установить новый поток вывода
    PYTHON_API void setStream(std::ostream &os) {
        std::lock_guard<std::mutex> lock(this->_mutex);
        this->_out = &os;
    }

    // Установить поток вывода в стандартный (stdout)
    PYTHON_API void resetStream() { this->setStream(std::cout); }

    // Добавить сообщение в буффер
    PYTHON_API void add(const std::string &msg) {
        std::lock_guard<std::mutex> lock(this->_mutex);
        this->_buffer.emplace_back("[" + now_timestamp_str() + "] " + msg);
    }

    // Сбросить буфер в поток вывода
    PYTHON_API void flush() {
        std::lock_guard<std::mutex> lock(this->_mutex);
        while (!_buffer.empty()) {
            (*_out) << _buffer.front() << '\n';
            _buffer.pop_front();
        }
        this->_out->flush();
    }
};

/**
 * RAII-Таймер для логирования времени выполнения блока кода.
 *
 * Используется в сочетании с Logger для записи времени начала и окончания выполнения блока кода.
 * Деструктор таймера автоматически посчитает время выполнения и запишет длительность исполнения.
 */
class PYTHON_API ScopedTimer {
  private:
    // имя таймера
    std::string _name;
    // ссылка на логгер
    Logger *_logger;
    // время запуска таймера
    std::chrono::high_resolution_clock::time_point _start;

  public:
    // Конструктор таймера с именем и ссылкой на логгер (по умолчанию синглтоновский), записывает старт
    PYTHON_API ScopedTimer(const std::string &name, Logger *logger = &Logger::instance())
        : _name(name), _logger(logger), _start(std::chrono::high_resolution_clock::now()) {
        if (_logger)
            _logger->add("[START] " + _name);
    }

    // Деструктор таймера, который посчитает время выполнения и запишет длительность исполнения.
    PYTHON_API ~ScopedTimer() {
        if (_logger) {
            auto end = std::chrono::high_resolution_clock::now();
            auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - _start).count();
            _logger->add("  [END] " + _name + " | duration: " + std::to_string(dur) + " μs");
        }
    }
};

// Отправить сообщение в логгер
#define LOG_MSG(msg) ::spindynapy::Logger::instance().add(msg);

// Отправить сообщение в логгер и сразу опубликовать
#define LOG_MSG_PRINT(msg)                                                                                   \
    ::spindynapy::Logger::instance().add(msg);                                                               \
    ::spindynapy::Logger::instance().flush();

// Определить таймер для логирования времени выполнения блока кода (деструктор посчитает время)
#define SCOPED_LOG_TIMER(name) ::spindynapy::ScopedTimer timer(name, &Logger::instance());

// Определить таймер и сразу печатать по его деструктуризации
#define SCOPED_LOG_TIMER_PRINT(name)                                                                         \
    ::spindynapy::ScopedTimer timer(name, &Logger::instance());                                              \
    ::spindynapy::Logger::instance().flush();

} // namespace spindynapy

// функция для связывания Логгера с Python
inline void pyBindLogger(py::module_ &module) {
    using namespace spindynapy;

    // АДАПТЕР для iostream, чтобы можно было передавать Python-объекты в Logger
    struct IOStreamAdapter {
        py::object keep_alive; // держим ссылку
        py::detail::pythonbuf buf;
        std::ostream os;
        explicit IOStreamAdapter(const py::object &o) : keep_alive(o), buf(keep_alive), os(&buf) {}
    };

    // Указатель-холдер. GC не удалит объект, nodelete.
    using NoDelete = std::unique_ptr<Logger, py::nodelete>;

    // -------- | LOGGER | --------
    py::module_ logger_module = module.def_submodule("logger");

    logger_module.doc() =
        "Логирование сообщений в stdout или в любой другой поток вывода.\n"
        "  Для отладки и мониторинга выполнения программы, длительности исполнения блоков кода.\n"
        "\n"
        "Используется для отладки и мониторинга выполнения программы.\n"
        "Вызывается через Logger::instance()";

    // холдер для сингтола, а то GC удалит объект
    static std::shared_ptr<IOStreamAdapter> singleton_logger_holder;

    py::class_<Logger, NoDelete>(logger_module, "Logger")
        .def(py::init([]() { return &Logger::instance(); }), py::return_value_policy::reference)
        .def_static(
            "instance",
            []() -> Logger & { return Logger::instance(); },
            py::return_value_policy::reference,
            py::doc("Получить объект логгера (синглтон)")
        )
        .def("add", &Logger::add, py::arg("msg"), py::doc("Добавить сообщение в буфер логгера"))
        .def("flush", &Logger::flush, py::doc("Сбросить буфер в поток вывода"))
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
            py::doc("Установить поток вывода для логгера")
        )

        .def(
            "reset_stream",
            [](Logger &self) {
                singleton_logger_holder.reset();
                self.resetStream();
            },
            py::doc("Сбросить поток вывода в стандартный (stdout)")
        )
        .doc() = "Логгер для вывода сообщений в поток.\n"
                 "\n"
                 "Используется для отладки и мониторинга выполнения программы.\n"
                 "Вызывается через Logger::instance()";

    py::class_<spindynapy::ScopedTimer>(logger_module, "ScopedTimer")
        .def(
            py::init<const std::string &, spindynapy::Logger *>(),
            py::arg("name"),
            py::arg("logger") = &Logger::instance()
        )
        .def(
            "__enter__", [](spindynapy::ScopedTimer &t) { return &t; }, py::doc("Запустить таймер")
        )
        .def(
            "__exit__", [](spindynapy::ScopedTimer &, py::args) {}, py::doc("Завершить таймер")
        )
        .doc() =
        "RAII-Таймер для логирования времени выполнения блока кода.\n"
        "\n"
        "Используется в сочетании с Logger для записи времени начала и окончания выполнения блока кода.\n"
        "Деструктор таймера автоматически посчитает время выполнения и запишет длительность исполнения.";
}

#endif // ! __LOGGER_HPP__
