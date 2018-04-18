#include <aspect/simulator.h>
#include <iostream>

/*
 * Launch the following function when this plugin is created. Use the update
 * script to update this parameter file. Unfortunately at this point the file
 * has already been loaded, so the current ASPECT instance would get the old
 * file. Thus we start a new ASPECT instance that will load the new file.
 * To avoid an endless recursion we remove the shared library from the new
 * input file (otherwise the new instance would call this library again, and so
 * on). After finishing the new instance we exit to not continue the old one.
 *
 * This test in particular calls the update script multiple times on the same
 * file, and ensures that only the first time changes happen. Every subsequent
 * application should not change the file any more.
 */
int f()
{
  // call ASPECT with "--" and pipe an existing input file into it.
  int ret;
  std::string command;

  command = ("cat update_script_2.x.prm | sed 's:set Additional shared libraries = ./libupdate_script_2.so::' | "
             "sed -f " ASPECT_SOURCE_DIR "/doc/update_prm_files_to_2.0.0.sed > output-update_script_2/updated.prm && "
             "cat output-update_script_2/updated.prm | sed -f " ASPECT_SOURCE_DIR "/doc/update_prm_files_to_2.0.0.sed > output-update_script_2/updated2.prm && "
             "diff output-update_script_2/updated.prm output-update_script_2/updated2.prm");
  std::cout << "Executing the update script:\n"
            << command
            << std::endl;
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  command = ("../aspect output-update_script_2/updated2.prm");
  std::cout << "Running ASPECT with updated parameter file:\n"
            << command
            << std::endl;
  ret = system (command.c_str());
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  // abort current process:
  exit (0);
  return 42;
}


// run this function by initializing a global variable by it
int i = f();
