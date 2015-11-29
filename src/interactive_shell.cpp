#include "interactive_shell.h"

#include <exception>
#include <iostream>

std::string InteractiveShell::PROMPT_ = ">=> ";

void InteractiveShell::register_command(std::string name,
                                        std::unique_ptr<ShellCommand> command)
{
  if (commands_.find(name) != commands_.end()) {
    throw std::logic_error("InteractiveShell::register_command: command \"" +
                           name + "\" already registered");
  } else {
    commands_.emplace(name, std::move(command));
  }
}

void InteractiveShell::run()
{
  while (true) {
    std::cout << PROMPT_;

    std::string line;
    std::cin >> line;

    if (line == "") continue;

    std::size_t pos_first_space = line.find(' ');
    std::string command;
    if (pos_first_space != std::string::npos) {
      command = line;

    } else {
      command = line.substr(0, pos_first_space);
    }

    std::vector<std::string> args;
    size_t current;
    size_t next = pos_first_space;
    do {
      current = next + 1;
      next = line.find(' ', current);
      args.emplace_back(line.substr(current, next - current));
    } while (next != std::string::npos);

    execute(command, args);
  }
}

void InteractiveShell::execute(std::string name,
                               const std::vector<std::string>& args)
{
  CommandMap::const_iterator cmd_pair = commands_.find(name);
  if (cmd_pair != commands_.end()) {
    std::cerr << "ERROR: Command \"" + name + "\" not registered\n";
  } else {
    cmd_pair->second->execute(args);
  }
}
