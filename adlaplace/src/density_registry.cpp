#include "adlaplace/api/density_registry.hpp"

#include <map>
#include <stdexcept>

static std::map<std::string, LogDensObsFn>& obs_registry() {
  static std::map<std::string, LogDensObsFn> reg;
  return reg;
}

static std::map<std::string, LogDensSingleDataFn>& single_data_registry() {
  static std::map<std::string, LogDensSingleDataFn> reg;
  return reg;
}

static std::map<std::string, LogDensSingleRandomDiagFn>& single_random_diag_registry() {
  static std::map<std::string, LogDensSingleRandomDiagFn> reg;
  return reg;
}

void register_log_dens_obs(const std::string& name, LogDensObsFn fn) {
  if (name.empty() || fn == nullptr) {
    throw std::runtime_error("register_log_dens_obs: invalid name or function pointer");
  }
  auto& reg = obs_registry();
  if (reg.count(name)) {
    throw std::runtime_error("register_log_dens_obs: name already registered: " + name);
  }
  reg[name] = fn;
}

void register_log_dens_single_data(const std::string& name, LogDensSingleDataFn fn) {
  if (name.empty() || fn == nullptr) {
    throw std::runtime_error("register_log_dens_single_data: invalid name or function pointer");
  }
  auto& reg = single_data_registry();
  if (reg.count(name)) {
    throw std::runtime_error("register_log_dens_single_data: name already registered: " + name);
  }
  reg[name] = fn;
}

void register_log_dens_single_random_diag(const std::string& name, LogDensSingleRandomDiagFn fn) {
  if (name.empty() || fn == nullptr) {
    throw std::runtime_error("register_log_dens_single_random_diag: invalid name or function pointer");
  }
  auto& reg = single_random_diag_registry();
  if (reg.count(name)) {
    throw std::runtime_error("register_log_dens_single_random_diag: name already registered: " + name);
  }
  reg[name] = fn;
}

LogDensObsFn resolve_log_dens_obs(const std::string& name) {
  const auto& reg = obs_registry();
  auto it = reg.find(name);
  if (it == reg.end()) {
    throw std::runtime_error("unknown log_dens_obs name: " + name);
  }
  return it->second;
}

LogDensSingleDataFn resolve_log_dens_single_data(const std::string& name) {
  const auto& reg = single_data_registry();
  auto it = reg.find(name);
  if (it == reg.end()) {
    throw std::runtime_error("unknown log_dens_single_data name: " + name);
  }
  return it->second;
}

LogDensSingleRandomDiagFn resolve_log_dens_single_random_diag(const std::string& name) {
  const auto& reg = single_random_diag_registry();
  auto it = reg.find(name);
  if (it == reg.end()) {
    throw std::runtime_error("unknown log_dens_single_random_diag name: " + name);
  }
  return it->second;
}

bool log_dens_single_uses_random_diag(const std::string& name) {
  return single_random_diag_registry().count(name) > 0;
}
