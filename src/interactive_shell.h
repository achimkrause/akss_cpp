#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

class ShellCommand
{
 public:
  virtual ~ShellCommand();
  virtual bool execute(const std::vector<std::string>& args) = 0;
};

using CommandMap = std::map<std::string, std::unique_ptr<ShellCommand>>;

class InteractiveShell
{
 public:
  void register_command(std::string name,
                        std::unique_ptr<ShellCommand> command);
  [[noreturn]] void run();

 private:
  void execute(std::string name, const std::vector<std::string>& args);

  static std::string PROMPT_;
  CommandMap commands_;
};
